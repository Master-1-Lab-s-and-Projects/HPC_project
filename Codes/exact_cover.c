#include "exact_cover.h"
#include "utility.h"
#include <stdio.h>
#include <stdlib.h>
#include <err.h>
#include <string.h>

#include <mpi.h>

#define WAS_CHOSEN 1
#define DISTRIBUTION_LEVEL_THRESHOLD 0

const int MPI_TAG_SOLUTION_COUNT = 1;   // tag pour la communication du nb de solutions (MPI)
const int MPI_TAG_STATUS_QUERY = 2;     // tag pour la demande d'état du context d'un proc (MPI)
const int MPI_TAG_DISTRIB_QUERY = 3;    // tag pour la demande de distribution du travail d'un proc (MPI)
const int MPI_TAG_NEW_WORK = 4;         // tag pour l'envoi d'un nouveau travail à un proc libre (MPI)


int *distribution_buffer = NULL;    //
int *status_buffer = NULL;          // contient la len, les child_num, et les upper_bounds du context
int **procs_status = NULL;          // contient les status de tous les procs esclaves

bool bypass = false;


void solution_found(const struct instance_t *instance, struct context_t *ctx)
{
//#pragma omp atomic
    ctx->solutions++;
    if (!print_solutions)
        return;
    printf("Trouvé une nouvelle solution au niveau %d après %lld noeuds\n",
            ctx->level, ctx->nodes);
    printf("Options : \n");
    for (int i = 0; i < ctx->level; i++) {
        int option = ctx->chosen_options[i];
        printf("+ %d : ", option);
        print_option(instance, option);
    }
    printf("\n");
    printf("----------------------------------------------------\n");
}

void cover(const struct instance_t *instance, struct context_t *ctx, int item);

void choose_option(const struct instance_t *instance, struct context_t *ctx,
        int option, int chosen_item, int was_chosen)
{
    ctx->chosen_options[ctx->level] = option;
    ctx->level++;
    for (int p = instance->ptr[option]; p < instance->ptr[option + 1]; p++) {
        int item = instance->options[p];
        if (item == chosen_item && was_chosen)
            continue;
        cover(instance, ctx, item);
    }
}

void uncover(const struct instance_t *instance, struct context_t *ctx, int item);

void unchoose_option(const struct instance_t *instance, struct context_t *ctx,
        int option, int chosen_item, int was_chosen)
{
    for (int p = instance->ptr[option + 1] - 1; p >= instance->ptr[option]; p--) {
        int item = instance->options[p];
        if (item == chosen_item && was_chosen)
            continue;
        uncover(instance, ctx, item);
    }
    ctx->level--;
}


/**
 * Returns the item associated with the active
 * option containing the least number of items
 * (the smallest active option).
 */
int choose_next_item(struct context_t *ctx)
{
    int best_item = -1;
    int best_options = 0x7fffffff;
    struct sparse_array_t *active_items = ctx->active_items;

    // TODO: small number of objects, check if to parallelize
    // active_items->len should be above ~=> 10 or so
    for (int i = 0; i < active_items->len; i++) {
        int item = active_items->p[i];
        struct sparse_array_t *active_options = ctx->active_options[item];
        int k = active_options->len;
        if (k < best_options) {
            best_item = item;
            best_options = k;
        }
    }
    return best_item;
}

void deactivate(const struct instance_t *instance, struct context_t *ctx,
        int option, int covered_item);

void cover(const struct instance_t *instance, struct context_t *ctx, int item)
{
    if (item_is_primary(instance, item))
        sparse_array_remove(ctx->active_items, item);
    struct sparse_array_t *active_options = ctx->active_options[item];
    for (int i = 0; i < active_options->len; i++) {
        int option = active_options->p[i];
        deactivate(instance, ctx, option, item);
    }
}


void deactivate(const struct instance_t *instance, struct context_t *ctx,
        int option, int covered_item)
{
    for (int k = instance->ptr[option]; k < instance->ptr[option+1]; k++) {
        int item = instance->options[k];
        if (item == covered_item)
            continue;
        sparse_array_remove(ctx->active_options[item], option);
    }
}


void reactivate(const struct instance_t *instance, struct context_t *ctx,
        int option, int uncovered_item);

void uncover(const struct instance_t *instance, struct context_t *ctx, int item)
{
    struct sparse_array_t *active_options = ctx->active_options[item];
    for (int i = active_options->len - 1; i >= 0; i--) {
        int option = active_options->p[i];
        reactivate(instance, ctx, option, item);
    }
    if (item_is_primary(instance, item))
        sparse_array_unremove(ctx->active_items);
}


void reactivate(const struct instance_t *instance, struct context_t *ctx,
        int option, int uncovered_item)
{
//#pragma omp parallel for
    for (int k = instance->ptr[option + 1] - 1; k >= instance->ptr[option]; k--) {
        int item = instance->options[k];
        if (item == uncovered_item)
            continue;
        sparse_array_unremove(ctx->active_options[item]);
    }
}

int master_queried(int tag)
{
    MPI_Status status;
    //MPI_Request request;

    //if (rank == 1 && tag == MPI_TAG_STATUS_QUERY && ctx->solutions > 1300000)
    //    printf("Proc [%d]: here\n", rank);

    bool query_value = false;
    int request_completed = 1;
    if (rank == 1) {
        if (MPI_Recv(&query_value, 1, MPI_C_BOOL, 0, tag, MPI_COMM_WORLD, &status)
                != MPI_SUCCESS)
            printf("Proc [%d]: query unsuccessful[_%d_]\n", rank, tag);

        //if (query_value)
        //    printf("Proc [%d]: query [_%d_]\n", rank, tag);

        //MPI_Test(&request, &request_completed, &status);
    }


    return request_completed && query_value;
}

/**
 * Sends upper_bounds, and child_num to the master
 */
void send_status(const struct instance_t *instance, const struct context_t *ctx)
{
    int n = instance->n_items;
    status_buffer[0] = ctx->level;
    memcpy(status_buffer + 1, ctx->child_num,
            (ctx->level + 1) * sizeof(*ctx->child_num));
    memcpy(status_buffer + 1 + n, ctx->upper_bounds,
            (ctx->level + 1) * sizeof(*ctx->upper_bounds));
    printf("Proc [%d]: sending STATUS to master\n", rank);

    MPI_Send(status_buffer, 2 * n + 1, MPI_INTEGER, 0,
            MPI_TAG_STATUS_QUERY, MPI_COMM_WORLD);
    printf("Proc [%d]: STATUS sent to master\n", rank);
}

int get_distribution_lvl(const struct context_t *ctx)
{
    int distrib_lvl = 0; // get the highest level not completed

    // TODO: fix ---> it has to be pretty
    if (ctx->child_num[0] + nb_proc - 1 >= ctx->upper_bounds[0])
        distrib_lvl++;

    while (ctx->child_num[distrib_lvl] + 1 >= ctx->upper_bounds[distrib_lvl])
        distrib_lvl++;

    return distrib_lvl;
}

void distribute_work(const struct instance_t *instance, struct context_t *ctx)
{
    int distrib_lvl = get_distribution_lvl(ctx);
    int distrib_l_bound = 0;

    if (distrib_lvl > DISTRIBUTION_LEVEL_THRESHOLD) {
        distrib_lvl = -1;
        bypass = true;
        printf("Proc [%d]: === WILL NOT DISTRIBUTE ANYMORE ===\n", rank);

    } else {
        distrib_l_bound = ctx->child_num[distrib_lvl]
            + ((distrib_lvl == 0) ? (nb_proc - 1) : 1);
        ctx->upper_bounds[distrib_lvl] = distrib_l_bound;
        memcpy(distribution_buffer + 2, ctx->chosen_options, distrib_lvl * sizeof(*distribution_buffer));

        printf("Proc [%d]: child_num[%d] = %d\n",
                rank, distrib_lvl, ctx->child_num[distrib_lvl]);
    }

    distribution_buffer[0] = distrib_lvl;
    distribution_buffer[1] = distrib_l_bound;

    printf("Proc [%d]: distributing level %d from child %d. New upper_bound: %d\n",
            rank, distrib_lvl, distrib_l_bound,
            (distrib_lvl == -1) ? -1 : ctx->upper_bounds[distrib_lvl]);

    MPI_Send(distribution_buffer, 2 + instance->n_items, MPI_INTEGER, 0,
            MPI_TAG_DISTRIB_QUERY, MPI_COMM_WORLD);
    printf("Proc [%d]: sent DISTRIB_INFO to master\n", rank);
}

/**
 * TODO: Use omp task + threadprivate (isolate context between threads)
 *
 * lower_bound and upper_bound indicate the bounds set for the current
 * level. If -1 is specified, the usual values (0 and active_options->len)
 * are used. It is useful for the level 0 (as the lower bound is linked to
 * the rank of the process) and when using the dynamic distribution,
 * to run the solver from a certain branch of the tree.
 *
 * Note: only the lower bound need to be specified.
 */
void solve(const struct instance_t * const instance,
        struct context_t * const ctx, int lower_bound)
{
    ctx->nodes++;
    if (ctx->nodes == next_report)
        progress_report(ctx);
    if (sparse_array_empty(ctx->active_items)) {
        solution_found(instance, ctx);
        return;                         /* succès : plus d'objet actif */
    }


    int chosen_item = choose_next_item(ctx);
    struct sparse_array_t *active_options = ctx->active_options[chosen_item];
    if (sparse_array_empty(active_options))
        return;           /* échec : impossible de couvrir chosen_item */

    cover(instance, ctx, chosen_item);
    ctx->num_children[ctx->level] = active_options->len;


    ctx->lower_bounds[ctx->level] = ((lower_bound == -1) ? 0 : lower_bound);
    ctx->upper_bounds[ctx->level] = active_options->len;

    if (ctx->child_num[0] + 2 * (nb_proc - 1) >= ctx->upper_bounds[0] && !bypass) {
        if (master_queried(MPI_TAG_STATUS_QUERY)) {
            printf("Proc [%d]: Master requested STATUS\n", rank);
            send_status(instance, ctx);
        }

        if (master_queried(MPI_TAG_DISTRIB_QUERY)) {
            printf("Proc [%d]: Master requested DISTRIBUTION\n", rank);
            distribute_work(instance, ctx);
        }
    }

    if (ctx->level == 0) {
        printf("Proc [%d]: l_bound = %d && u_bound = %d\n",
                rank, ctx->lower_bounds[ctx->level], ctx->upper_bounds[ctx->level]);
    }

    for (int k = ctx->lower_bounds[ctx->level]; k < ctx->upper_bounds[ctx->level];) {
        //if (ctx->level == 0 && rank == 1)
        //    printf("Proc [%d]: for looop, k = %d\n", rank, k);
        int option = active_options->p[k];
        ctx->child_num[ctx->level] = k;
        choose_option(instance, ctx, option, chosen_item, WAS_CHOSEN);
        // On assigne les taches ici, -> idée savoir qu'il ont fini leur travail quand il sort de cette boucle avec l'id k de la ligne 564
        // /!\ On envoie au tache le contexte, après avoir choisi l'option, et on sort de la boucle quand on a fini solve au niveau k :peepoez:
        // faudra faire marcher ca

        //#pragma omp task
        solve(instance, ctx, -1);

        if (ctx->solutions >= max_solutions)
            return;
        //apres on communique le nb de solutions et on attends une prochaine tacches
        // Si on veut aller plus loin -> methode arbre binomial comme ca on parallelise le plus possible mais il risque d'avoir plus de noeuds travailleurs
        unchoose_option(instance, ctx, option, chosen_item, WAS_CHOSEN);


        // On ne distribue qu'au premier niveau
        k += (ctx->level == 0) ? (nb_proc - 1) : 1;
    }
    if (ctx->level == 0) {
        printf("Proc [%d]: l_bound = %d && u_bound = %d\n",
                rank, ctx->lower_bounds[ctx->level], ctx->upper_bounds[ctx->level]);
    }

    //Typiquement on rajoute une boucle ici pour les proc 1 ... n-1 et on les fait tourner a l'infini, on les arrete avec un signal depuis le programme principal

    uncover(instance, ctx, chosen_item);                      /* backtrack */
}



int wait_for_done_proc(long long *rcvd_solutions)
{
    MPI_Status status;
    MPI_Recv(rcvd_solutions, 1, MPI_LONG_LONG_INT, MPI_ANY_SOURCE,
            MPI_TAG_SOLUTION_COUNT, MPI_COMM_WORLD, &status);
    return status.MPI_SOURCE;
}

int choose_proc_to_distribute(int proc_to_exclude, int nb_items)
{
    int lowest_common_lvl = nb_items - 1;
    for (int i=1; i < nb_proc; i++) {
        if (i == proc_to_exclude)
            continue;

        memset(procs_status[i - 1], 0, (2 * nb_items + 1) * sizeof(*procs_status[0]));

        printf("Master: query STATUS of proc [%d]\n", i);
        // Collect the status of each process.
        bool query_status = true;
        MPI_Send(&query_status, 1, MPI_C_BOOL,
                i, MPI_TAG_STATUS_QUERY, MPI_COMM_WORLD);

        if (MPI_Recv(procs_status[i - 1], 2 * nb_items + 1, MPI_INTEGER, i,
                    MPI_TAG_STATUS_QUERY, MPI_COMM_WORLD, MPI_STATUS_IGNORE)
                != MPI_SUCCESS)
            printf("Master: reception error");

        printf("Master: recvd status of proc [%d]\n", i);
        // Determine the lowest level in the
        // backtracking tree reached by all the process.
        int current_lvl = procs_status[i - 1][0];
        printf("\tCurrent level of proc [%d]: %d\n", i, current_lvl);
        if (current_lvl <= lowest_common_lvl) {
            lowest_common_lvl = current_lvl;
        }
    }


    int highest_not_completed_lvl = lowest_common_lvl;
    int proc_to_distribute = -1;
    int max_children_to_expl = 0;
    for (int i=1; i < nb_proc; i++) {
        if (i == proc_to_exclude)
            continue;

        int *child_num = &procs_status[i - 1][1];
        int *upper_bounds = &procs_status[i - 1][1 + nb_items];
        int lvl = 0;

        // Find the highest not completed level for the proc of rank i.
        // Meaning, the process is currently exploring the last configuration
        // of options for all the levels above.
        while (lvl < highest_not_completed_lvl && child_num[lvl] == upper_bounds[lvl] - 1)
            lvl++;


        int children_to_expl = upper_bounds[lvl] - child_num[lvl] - 1;
        if (lvl < highest_not_completed_lvl
                || (lvl == highest_not_completed_lvl
                    && children_to_expl > max_children_to_expl)) {

            highest_not_completed_lvl = lvl;
            proc_to_distribute = i;
            max_children_to_expl = children_to_expl;
        }
    }

    /** TODO add a threshold. Example
     *  if (highest_not_completed_lvl > A_CER TAIN_DEPTH) // based on the number of items
     *      proc_to_distribute = -1; // equiv. to: do not distribute
     */

    return proc_to_distribute;
}

void distribute_work_master(const struct instance_t *instance,
        int worker, int *distribution_buffer)
{
    /**
     * To distribute work we need to get:
     * - the level 'l' to distribute
     * - the lower bound for the level 'l'
     * - the first 'l' chosen options
     */

    bool query_value = true;
    MPI_Send(&query_value, 1, MPI_C_BOOL, worker,
            MPI_TAG_DISTRIB_QUERY, MPI_COMM_WORLD);
    printf("Master: sent DISTRIB_QUERY to proc [%d]\n", worker);

    MPI_Recv(distribution_buffer, 2 + instance->n_items, MPI_INTEGER, worker,
            MPI_TAG_DISTRIB_QUERY, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    printf("Master: rcvd DISTRIB_QUERY response from proc [%d]\n", worker);

}

void build_context(const struct instance_t *instance,
        struct context_t *ctx, const int *chosen_options, int n_options)
{
    ctx->solutions = 0;
    // We don't reset the number of nodes browsed. This allows us
    // to check that all the branches have been explored.

    for (int i=0; i < n_options; i++) {
        choose_option(instance, ctx, chosen_options[i], -1, !WAS_CHOSEN);
    }
    printf("Proc [%d]: Context built, level = %d\n", rank, ctx->level);
}

void reset_context(const struct instance_t *instance,
        struct context_t *ctx, const int *chosen_options, int n_options)
{
    for (int i=n_options - 1; i >= 0; i--) {
        unchoose_option(instance, ctx, chosen_options[i], -1, !WAS_CHOSEN);
    }
}

void launch_parallel(const struct instance_t *instance, struct context_t *ctx)
{
    if (rank == 0) {
        printf("Number of proc: %d\n", nb_proc);

        int nb_free_worker = 0;
        bool distribute_remaining_work = true;
        while (nb_free_worker < nb_proc - 1) {
            long long rcvd_solutions;
            int done_proc;

            // collect the number of solutions
            done_proc = wait_for_done_proc(&rcvd_solutions);
            ctx->solutions += rcvd_solutions;
            nb_free_worker++;

            printf("Master: Received %lld solutions from proc %d. "
                    "Total solutions: %lld.\n",
                    rcvd_solutions, done_proc, ctx->solutions);

            if (distribute_remaining_work) {
                int proc_to_distribute;
                proc_to_distribute = choose_proc_to_distribute(done_proc, instance->n_items);
                printf("Master: proc n°%d will distribute its work\n", proc_to_distribute);

                distribute_work_master(instance, proc_to_distribute, distribution_buffer);

                if (distribution_buffer[0] == -1) {
                    distribute_remaining_work = false;
                    printf("Master: === END of DISTRIBUTION ===\n");
                } else {
                    nb_free_worker--;
                }
            }

            MPI_Send(distribution_buffer, 2 + instance->n_items, MPI_INTEGER,
                    done_proc, MPI_TAG_NEW_WORK, MPI_COMM_WORLD);
            printf("Master: sent NEW_WORK to proc [%d]\n", done_proc);
        }

        printf("FINI. Trouvé %lld solutions en %.1fs\n", ctx->solutions,
                wtime() - start);

    } else {
        int lower_bound = rank - 1;
        int n_options = -1;
        distribution_buffer[0] = -1;

        // loop to re-launch when done
        do {
            if (distribution_buffer[0] != -1) {
                n_options = distribution_buffer[0];
                lower_bound = distribution_buffer[1];
                build_context(instance, ctx, distribution_buffer + 2, n_options);
            }

            printf("Proc [%d]: [%lld] solutions before solve\n",
                    rank, ctx->solutions);
            solve(instance, ctx, lower_bound);
            printf("Proc [%d]: done with [%lld] solutions\n",
                    rank, ctx->solutions);

            MPI_Send(&ctx->solutions, 1, MPI_LONG_LONG_INT, 0,
                    MPI_TAG_SOLUTION_COUNT, MPI_COMM_WORLD);

            if (distribution_buffer[0] != -1) {
                reset_context(instance, ctx, distribution_buffer + 2, n_options);
            }

            MPI_Recv(distribution_buffer, 2 + instance->n_items, MPI_INTEGER, 0,
                    MPI_TAG_NEW_WORK, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            printf("Proc [%d], next depth [%d], lower_bound [%d]\n",
                    rank, distribution_buffer[0], distribution_buffer[1]);

        } while(distribution_buffer[0] != -1);

        printf("Proc [%d]: DONE\n", rank);
    }
}


/**
 * TODO fix error messages
 */
void parallel_setup(const struct instance_t *instance)
{
    int n = instance->n_items;
    if (rank == 0) {
        procs_status = malloc((nb_proc - 1) * sizeof(*procs_status));
        if (procs_status == NULL)
            err(1, "impossible d'allouer le buffer de status");

        for (int i=0; i < nb_proc - 1; i++) {
            procs_status[i] = malloc((2 * n + 1) * sizeof(*procs_status[i]));
            if (procs_status[i] == NULL)
                err(1, "impossible d'allouer le buffer de status");
        }

    } else {
        status_buffer = malloc((2 * n + 1) * sizeof(*status_buffer));
        if (status_buffer == NULL)
            err(1, "impossible d'allouer le buffer de status");
    }

    distribution_buffer = malloc((n + 2) * sizeof(*distribution_buffer));
    if (distribution_buffer == NULL)
        err(1, "impossible d'allouer le buffer des options choisies");
}

struct context_t * backtracking_setup(const struct instance_t *instance)
{
    struct context_t *ctx = malloc(sizeof(*ctx));
    if (ctx == NULL)
        err(1, "impossible d'allouer un contexte");
    ctx->level = 0;
    ctx->nodes = 0;
    ctx->solutions = 0;
    int n = instance->n_items;
    int m = instance->n_options;
    ctx->active_options = malloc(n * sizeof(*ctx->active_options));
    ctx->chosen_options = malloc(n * sizeof(*ctx->chosen_options));
    ctx->child_num = malloc(n * sizeof(*ctx->child_num));
    ctx->num_children = malloc(n * sizeof(*ctx->num_children));
    ctx->lower_bounds = malloc(n * sizeof(*ctx->lower_bounds));
    ctx->upper_bounds = malloc(n * sizeof(*ctx->upper_bounds));

    status_buffer = malloc((2 * n + 1) * sizeof(*status_buffer));

    if (ctx->active_options == NULL || ctx->chosen_options == NULL
            || ctx->child_num == NULL || ctx->num_children == NULL
            || ctx->lower_bounds == NULL || ctx->upper_bounds == NULL
            || status_buffer == NULL)
        err(1, "impossible d'allouer le contexte");

    memset(ctx->child_num, 0, n * sizeof(*ctx->child_num));
    memset(ctx->num_children, 0, n * sizeof(*ctx->num_children));
    memset(ctx->lower_bounds, 0, n * sizeof(*ctx->lower_bounds));
    memset(ctx->upper_bounds, 0, n * sizeof(*ctx->upper_bounds));


    ctx->active_items = sparse_array_init(n);

    // TODO: probably not possible... modification to array
    for (int item = 0; item < instance->n_primary; item++)
        sparse_array_add(ctx->active_items, item);

    // TODO: probably not possible... parallelized malloc !?
    for (int item = 0; item < n; item++)
        ctx->active_options[item] = sparse_array_init(m);

    /*
     * [o_1, ..., 0_m]
     * instance->ptr[o_i] <= k < instance->ptr[o_(i + 1)]
     * TODO: check options disjointed..
     */
    for (int option = 0; option < m; option++)
        for (int k = instance->ptr[option]; k < instance->ptr[option + 1]; k++) {
            int item = instance->options[k];
            sparse_array_add(ctx->active_options[item], option);
        }

    return ctx;
}
