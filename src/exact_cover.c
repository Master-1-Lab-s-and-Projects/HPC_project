#include "exact_cover.h"
#include "utility.h"
#include <stdio.h>
#include <stdlib.h>
#include <err.h>
#include <string.h>
#include <assert.h>

#include <mpi.h>

#ifdef DEBUG
    #define DPRINTF(...) printf(__VA_ARGS__)
#else
    #define DPRINTF(...)
#endif

#define min(a,b) \
    ({ __typeof__ (a) _a = (a); \
     __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

#define distribution_level *distribution_info.level
#define distribution_lower_bound *distribution_info.lower_bound
#define distribution_chosen_options distribution_info.chosen_options

#define WAS_CHOSEN 1
#define DISTRIBUTION_LEVEL_THRESHOLD 1

const int MPI_TAG_STATUS_QUERY = 0;     // tag pour la demande d'état du context d'un proc (MPI)
const int MPI_TAG_DISTRIB_QUERY = 1;    // tag pour la demande de distribution du travail d'un proc (MPI)
const int MPI_TAG_SOLUTION_COUNT = 2;   // tag pour la communication du nb de solutions (MPI)
const int MPI_TAG_NEW_WORK = 3;         // tag pour l'envoi d'un nouveau travail à un proc libre (MPI)


int *distribution_buffer = NULL;    //
int *status_buffer = NULL;          // contient la len, les child_num, et les upper_bounds du context
int **procs_status = NULL;          // contient les status de tous les procs esclaves

MPI_Request *state_reqs = NULL; // maybe state_requests
int *cmpl_state_reqs_indices = NULL;

MPI_Request master_reqs[2];
bool Irecv_called[2] = {false, false};
bool query_value = false;

bool stop_distribution = false;

struct distribution_info {
    int *level;
    int *lower_bound;
    int *chosen_options;
};
struct distribution_info distribution_info = {0};


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
        if (k <= best_options) {
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

int master_requested(int tag)
{

    //if (rank == 1 && tag == MPI_TAG_STATUS_QUERY)
    //    DPRINTF("Proc [%d]: here\n", rank);

    if (!Irecv_called[tag]) {
        MPI_Irecv(&query_value, 1, MPI_C_BOOL, 0, tag, MPI_COMM_WORLD, &master_reqs[tag]);
        Irecv_called[tag] = true;
    }
    //MPI_Recv(&query_value, 1, MPI_C_BOOL, 0, tag, MPI_COMM_WORLD, &status);

    int request_completed = 0;
    MPI_Test(&master_reqs[tag], &request_completed, MPI_STATUS_IGNORE);

    if (request_completed) {
        Irecv_called[tag] = false;
    }

    return request_completed;
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
    DPRINTF("Proc [%d]: sending STATUS to master\n", rank);

    MPI_Send(status_buffer, 2 * n + 1, MPI_INTEGER, 0,
            MPI_TAG_STATUS_QUERY, MPI_COMM_WORLD);
    DPRINTF("Proc [%d]: STATUS sent to master\n", rank);
}

/**
 * TODO: clean up
 */
// Find the highest not completed level for the proc of rank i.
// Meaning, the process is currently exploring the last configuration
// of options for all the levels above.
int get_distribution_lvl(int *child_num, int* upper_bound, int lvl_limit)
{
    int distrib_lvl = 0; // get the highest level not completed

    if (distrib_lvl < lvl_limit
            && child_num[0] + nb_proc - 1 >= upper_bound[0])
        distrib_lvl++;

    while (distrib_lvl < lvl_limit
            && child_num[distrib_lvl] + 1 >= upper_bound[distrib_lvl])
        distrib_lvl++;

    return distrib_lvl;
}


void distribute_work(const struct instance_t *instance, struct context_t *ctx)
{
    int distrib_lvl = get_distribution_lvl( ctx->child_num, ctx->upper_bounds, instance->n_items);

    if (distrib_lvl > DISTRIBUTION_LEVEL_THRESHOLD) {
        DPRINTF("Proc [%d]: === WILL NOT DISTRIBUTE ANYMORE (%d > %d) ===\n",
                rank, distrib_lvl, DISTRIBUTION_LEVEL_THRESHOLD);
        distrib_lvl = -1;
        stop_distribution = true;

    } else {
        distribution_lower_bound = ctx->child_num[distrib_lvl] + ((distrib_lvl == 0) ? (nb_proc - 1) : 1);

        ctx->upper_bounds[distrib_lvl] = distribution_lower_bound;
        memcpy(distribution_chosen_options, ctx->chosen_options,
                distrib_lvl * sizeof(*ctx->chosen_options));

        DPRINTF("Proc [%d]: child_num[%d] = %d\n",
                rank, distrib_lvl, ctx->child_num[distrib_lvl]);
    }

    distribution_level = distrib_lvl;

    DPRINTF("Proc [%d]: distributing level %d from child %d. New upper_bound: %d\n",
            rank, distrib_lvl, distribution_lower_bound,
            (distrib_lvl == -1) ? -1 : ctx->upper_bounds[distrib_lvl]);

    MPI_Send(distribution_buffer, 2 + instance->n_items, MPI_INTEGER, 0,
            MPI_TAG_DISTRIB_QUERY, MPI_COMM_WORLD);
    DPRINTF("Proc [%d]: sent DISTRIB_INFO to master\n", rank);
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
    int chosen_item_before = chosen_item;
    struct sparse_array_t *active_options = ctx->active_options[chosen_item];
    if (sparse_array_empty(active_options))
        return;           /* échec : impossible de couvrir chosen_item */
    //if (ctx->level == 0)
    //    DPRINTF("Proc [%d]: lvl[%d], chosen item = %d\n",
    //            rank, ctx->level, chosen_item);

    cover(instance, ctx, chosen_item);
    ctx->num_children[ctx->level] = active_options->len;
    ctx->lower_bounds[ctx->level] = ((lower_bound == -1) ? 0 : lower_bound);
    ctx->upper_bounds[ctx->level] = active_options->len;

    if (!stop_distribution && ctx->nodes % 100000 == 0) {
        if (master_requested(MPI_TAG_STATUS_QUERY)) {
            DPRINTF("Proc [%d]: Master requested STATUS\n", rank);
            send_status(instance, ctx);
        }

        if (master_requested(MPI_TAG_DISTRIB_QUERY)) {
            DPRINTF("Proc [%d]: Master requested DISTRIBUTION\n", rank);
            distribute_work(instance, ctx);
        }
    }

    //if (ctx->level == 0) {
    //    DPRINTF("Proc [%d]: l_bound = %d && u_bound = %d\n",
    //            rank, ctx->lower_bounds[ctx->level], ctx->upper_bounds[ctx->level]);
    //}

    for (int k = ctx->lower_bounds[ctx->level]; k < ctx->upper_bounds[ctx->level];) {
        //if (ctx->level == 0)
        //    DPRINTF("Proc [%d]: for looop k = %d\n", rank, k);
        int option = active_options->p[k];
        ctx->child_num[ctx->level] = k;
        choose_option(instance, ctx, option, chosen_item, WAS_CHOSEN);

        solve(instance, ctx, -1);

        if (ctx->solutions >= max_solutions)
            return;

        unchoose_option(instance, ctx, option, chosen_item, WAS_CHOSEN);


        // On ne distribue qu'au premier niveau
        k += (ctx->level == 0) ? (nb_proc - 1) : 1;
    }

    //if (ctx->level == 0) {
    //    DPRINTF("Proc [%d]: l_bound = %d && u_bound = %d\n",
    //            rank, ctx->lower_bounds[ctx->level], ctx->upper_bounds[ctx->level]);
    //}

    //Typiquement on rajoute une boucle ici pour les proc 1 ... n-1 et on les fait tourner a l'infini, on les arrete avec un signal depuis le programme principal

    uncover(instance, ctx, chosen_item);                      /* backtrack */
    assert(choose_next_item(ctx) == chosen_item_before);
}



int wait_for_done_proc(long long *rcvd_solutions)
{
    MPI_Status status;
    MPI_Recv(rcvd_solutions, 1, MPI_LONG_LONG_INT, MPI_ANY_SOURCE,
            MPI_TAG_SOLUTION_COUNT, MPI_COMM_WORLD, &status);
    return status.MPI_SOURCE;
}

/**
 * IDEA: Send to all workers the current time.
 * When a worker tests if a request has been received, if too much
 * time elapsed since the request has been issued (the time sent),
 * the request is ignored.
 */
void request_workers_status(int free_worker, int nb_items)
{
    int nb_reqs = 0;
    for (int i=1; i < nb_proc; i++) {
        if (i == free_worker)
            continue;

        procs_status[i][0] = -1;
        DPRINTF("Master: query STATUS of proc [%d]\n", i);

        //MPI_Request send_req;
        bool query_status = true;
        MPI_Send(&query_status, 1, MPI_C_BOOL, i,
                MPI_TAG_STATUS_QUERY, MPI_COMM_WORLD);

        MPI_Irecv(procs_status[i], 1 + 2 * nb_items, MPI_INTEGER, i,
                MPI_TAG_STATUS_QUERY, MPI_COMM_WORLD, &state_reqs[nb_reqs]);

        nb_reqs++;
    }

    DPRINTF("Master: waiting for workers to send state\n");
    int idx_count = 0;
    double curr_time = wtime();
    // TODO: add a define for the TIMEOUT
    while (wtime() - curr_time < 1 && idx_count == 0) {
        MPI_Testsome(nb_reqs, state_reqs, &idx_count,
                cmpl_state_reqs_indices, MPI_STATUSES_IGNORE);
    }

    DPRINTF("Master: %d processes returned\n", idx_count);
}

int choose_proc_to_distribute(int free_worker, int nb_items)
{
    int highest_not_completed_lvl = nb_items - 1;
    int proc_to_distribute = -1;
    int max_children_to_explore = 0;
    for (int i=1; i < nb_proc; i++) {
        int proc_current_lvl = procs_status[i][0];
        if (i == free_worker || proc_current_lvl == -1)
            continue;

        int *child_num = &procs_status[i][1];
        int *upper_bound = &procs_status[i][1 + nb_items];

        int lvl = get_distribution_lvl(child_num, upper_bound,
                min(proc_current_lvl, highest_not_completed_lvl));

        int children_to_explore = upper_bound[lvl] - child_num[lvl] - 1;
        if (lvl <= highest_not_completed_lvl
                && children_to_explore >= max_children_to_explore) {

            highest_not_completed_lvl = lvl;
            proc_to_distribute = i;
            max_children_to_explore = children_to_explore;
        }
    }

    /** TODO add a threshold. Example
     *  if (highest_not_completed_lvl > A_CER TAIN_DEPTH) // based on the number of items
     *      proc_to_distribute = -1; // equiv. to: do not distribute
     */

    return proc_to_distribute;
}

void request_work_distribution(const struct instance_t *instance,
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
    DPRINTF("Master: sent DISTRIB_QUERY to proc [%d]\n", worker);


    MPI_Request request;
    MPI_Irecv(distribution_buffer, 2 + instance->n_items, MPI_INTEGER, worker,
            MPI_TAG_DISTRIB_QUERY, MPI_COMM_WORLD, &request);

    double curr_time = wtime();
    int completed = 0;
    while (wtime() - curr_time < 2 && !completed) {
        MPI_Test(&request, &completed, MPI_STATUS_IGNORE);
    }
    if (completed)
        DPRINTF("Master: rcvd DISTRIB_QUERY response from proc [%d]\n", worker);
    else
        DPRINTF("Master: proc [%d] did not send DISTRIB_QUERY response\n", worker);

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
    DPRINTF("Proc [%d]: Context built, level = %d\n", rank, ctx->level);
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
        DPRINTF("Number of proc: %d\n", nb_proc);
        int nb_free_worker = 0;

        while (nb_free_worker < nb_proc - 1) {
            long long rcvd_solutions;
            int free_worker;

            // collect the number of solutions
            free_worker = wait_for_done_proc(&rcvd_solutions);
            ctx->solutions += rcvd_solutions;
            nb_free_worker++;

            DPRINTF("Master: Received %lld solutions from proc %d. "
                    "Total solutions: %lld.\n",
                    rcvd_solutions, free_worker, ctx->solutions);

            distribution_level = -1;
            if (nb_free_worker == 1) {
                int worker_to_distribute;
                request_workers_status(free_worker, instance->n_items);
                worker_to_distribute = choose_proc_to_distribute(free_worker, instance->n_items);
                DPRINTF("Master: proc n°%d will distribute its work\n", worker_to_distribute);

                if (worker_to_distribute == -1) {
                    DPRINTF("Master: === END of DISTRIBUTION ===\n");
                } else {
                    request_work_distribution(instance, worker_to_distribute, distribution_buffer);
                    if (distribution_level != -1)
                        nb_free_worker--;
                }
            }


            MPI_Send(distribution_buffer, 2 + instance->n_items, MPI_INTEGER,
                    free_worker, MPI_TAG_NEW_WORK, MPI_COMM_WORLD);
            DPRINTF("Master: sent NEW_WORK to proc [%d], free workers: %d\n",
                    free_worker, nb_free_worker);
        }

        printf("FINI. Trouvé %lld solutions en %.1fs\n", ctx->solutions,
                wtime() - start);

    } else {
        int lower_bound = rank - 1;
        distribution_level = -1;

        // loop to re-launch when done
        do {
            stop_distribution = false;
            if (distribution_level != -1) {
                lower_bound = distribution_lower_bound;
                build_context(instance, ctx, distribution_chosen_options, distribution_level);
            }

            int distrib_lvl_sav = distribution_level;
            solve(instance, ctx, lower_bound);
            distribution_level = distrib_lvl_sav;
            DPRINTF("Proc [%d]: done with [%lld] solutions\n", rank, ctx->solutions);

            MPI_Send(&ctx->solutions, 1, MPI_LONG_LONG_INT, 0, MPI_TAG_SOLUTION_COUNT, MPI_COMM_WORLD);

            if (distribution_level != -1) {
                reset_context(instance, ctx, distribution_chosen_options, distribution_level);
            }

            MPI_Recv(distribution_buffer, 2 + instance->n_items, MPI_INTEGER, 0,
                    MPI_TAG_NEW_WORK, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            DPRINTF("Proc [%d], next level [%d], lower_bound [%d]\n",
                    rank, distribution_level, distribution_lower_bound);

        } while(distribution_level != -1);

        DPRINTF("Proc [%d]: DONE\n", rank);
    }
}


/**
 * TODO fix error messages
 */
void parallel_setup(const struct instance_t *instance)
{
    int n = instance->n_items;
    if (rank == 0) {
        procs_status = malloc((nb_proc) * sizeof(*procs_status));
        if (procs_status == NULL)
            err(1, "impossible d'allouer le buffer de status");

        for (int i=0; i < nb_proc; i++) {
            procs_status[i] = malloc((1 + 2 * n) * sizeof(*procs_status[i]));
            if (procs_status[i] == NULL)
                err(1, "impossible d'allouer le buffer de status");
        }

        state_reqs = malloc(nb_proc * sizeof(*state_reqs));
        cmpl_state_reqs_indices = malloc(nb_proc * sizeof(*cmpl_state_reqs_indices));
        if (state_reqs == NULL || cmpl_state_reqs_indices == NULL)
            err(1, "impossible d'allouer les tableaux de req");


    } else {
        status_buffer = malloc((2 * n + 1) * sizeof(*status_buffer));
        if (status_buffer == NULL)
            err(1, "impossible d'allouer le buffer de status");
    }

    distribution_buffer = malloc((n + 2) * sizeof(*distribution_buffer));
    if (distribution_buffer == NULL)
        err(1, "impossible d'allouer le buffer des options choisies");

    distribution_info.level = &distribution_buffer[0];
    distribution_info.lower_bound = &distribution_buffer[1];
    distribution_info.chosen_options = &distribution_buffer[2];
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
