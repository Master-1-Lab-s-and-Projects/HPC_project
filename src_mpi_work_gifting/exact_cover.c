#include "exact_cover.h"
#include "utility.h"
#include <stdio.h>
#include <stdlib.h>
#include <err.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include <mpi.h>

#define INITIAL_NB_WORK_ORDERS 2 * nb_proc
const long long WORK_SHARE_DELTA = 1e7; // partage son travail tous les ... noeuds
const int MAX_LEVEL_SHARED = 50;


void solution_found(const struct instance_t *instance, struct context_t *ctx)
{
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
        int option, int chosen_item)
{
    ctx->chosen_options[ctx->level] = option;
    ctx->level++;
    for (int p = instance->ptr[option]; p < instance->ptr[option + 1]; p++) {
        int item = instance->options[p];
        if (item == chosen_item)
            continue;
        cover(instance, ctx, item);
    }
}

void uncover(const struct instance_t *instance, struct context_t *ctx, int item);

void unchoose_option(const struct instance_t *instance, struct context_t *ctx,
        int option, int chosen_item)
{
    for (int p = instance->ptr[option + 1] - 1; p >= instance->ptr[option]; p--) {
        int item = instance->options[p];
        if (item == chosen_item)
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
    for (int k = instance->ptr[option + 1] - 1; k >= instance->ptr[option]; k--) {
        int item = instance->options[k];
        if (item == uncovered_item)
            continue;
        sparse_array_unremove(ctx->active_options[item]);
    }
}

/**
 * Find the highest not completed level for the proc of rank i.
 * Meaning, the process is currently exploring the last configuration
 * of options for all the levels above.
 */
int get_distribution_lvl(int *child_num, int* num_children, int lvl_limit)
{
    int distrib_lvl = 0; // the highest level not completed

    if (distrib_lvl < lvl_limit && child_num[0] + nb_proc - 1 >= num_children[0])
        distrib_lvl++;

    while (distrib_lvl < lvl_limit && child_num[distrib_lvl] + 1 >= num_children[distrib_lvl])
        distrib_lvl++;

    return (distrib_lvl <= MAX_LEVEL_SHARED) ? distrib_lvl : -1;
}

int get_first_child_to_share(const struct context_t *ctx, int distrib_lvl)
{
    int remaining_children = ctx->num_children[distrib_lvl] - ctx->child_num[distrib_lvl];
    assert(remaining_children != 0);
    if (distrib_lvl == 0) {
        remaining_children -= remaining_children % (nb_proc - 1);
        remaining_children /= nb_proc - 1;
    }

    int new_remaining_children = ceil(((double) remaining_children)/ 2.0);

    if (distrib_lvl == 0)
        new_remaining_children *= nb_proc - 1;

    return ctx->child_num[distrib_lvl] + new_remaining_children;
}

void share_work(struct context_t *ctx, int *work_order_buffer, int work_order_size)
{
    ctx->next_work_share += WORK_SHARE_DELTA;

    int distrib_lvl = get_distribution_lvl(ctx->child_num, ctx->num_children, ctx->level);
    if (distrib_lvl == -1) {
        DPRINTF("Proc [%d]: Did not share work\n", rank);
        return;
    }

    int first_child_shared = get_first_child_to_share(ctx, distrib_lvl);
    int last_child_shared = ctx->num_children[distrib_lvl];
    DPRINTF("Proc [%d]: at level %d, [%d, ..., %d[ --> [%d, ..., %d[\n",
            rank, distrib_lvl,
            ctx->child_num[distrib_lvl], ctx->num_children[distrib_lvl],
            ctx->child_num[distrib_lvl], first_child_shared);

    memcpy(work_order_buffer + 2, ctx->child_num, (distrib_lvl + 1) * sizeof(*ctx->child_num));
    work_order_buffer[0] = distrib_lvl;
    work_order_buffer[1] = last_child_shared;
    work_order_buffer[2 + distrib_lvl] = first_child_shared;
    MPI_Send(work_order_buffer, work_order_size, MPI_INTEGER, 0, 0, MPI_COMM_WORLD);

    ctx->num_children[distrib_lvl] = first_child_shared;
}

void solve(const struct instance_t *instance, struct context_t *ctx,
        int *work_order_buffer, int work_order_size, int loop_incr,
        int distrib_lvl, int distrib_num_children)
{
    ctx->nodes++;
    if (ctx->nodes == next_report)
        progress_report(ctx);
    if (ctx->nodes == ctx->next_work_share)
        share_work(ctx, work_order_buffer, work_order_size);
    if (sparse_array_empty(ctx->active_items)) {
        solution_found(instance, ctx);
        return;                         /* succès : plus d'objet actif */
    }
    int chosen_item = choose_next_item(ctx);
    struct sparse_array_t *active_options = ctx->active_options[chosen_item];
    if (sparse_array_empty(active_options))
        return;           /* échec : impossible de couvrir chosen_item */
    cover(instance, ctx, chosen_item);

    if (ctx->level == distrib_lvl && distrib_num_children != -1)
        ctx->num_children[ctx->level] = distrib_num_children;
    else
        ctx->num_children[ctx->level] = active_options->len;

    for (int k = ctx->first_child[ctx->level]; k < ctx->num_children[ctx->level]; k += loop_incr) {
        if (ctx->level == 0 && 0)
            DPRINTF("#%d Proc [%d]\n", k, rank);
        int option = active_options->p[k];
        ctx->child_num[ctx->level] = k;
        choose_option(instance, ctx, option, chosen_item);
        solve(instance, ctx, work_order_buffer, work_order_size, 1,
                distrib_lvl, distrib_num_children);
        if (ctx->solutions >= max_solutions)
            return;
        unchoose_option(instance, ctx, option, chosen_item);
    }
    /* Reset the first child to be explored for the current level.  */
    ctx->first_child[ctx->level] = 0;

    uncover(instance, ctx, chosen_item);                      /* backtrack */
}

int ** allocate_work_orders(int **work_orders, int old_size, int new_size,
        int work_order_size)
{
    assert(old_size < new_size);
    work_orders = realloc(work_orders, new_size * sizeof(*work_orders));
    if (work_orders == NULL)
        err(1, "impossible d'allouer le tableau des ordres de travail");

    for (int i = old_size; i < new_size; i++) {
        work_orders[i] = malloc(work_order_size * sizeof(work_orders[i]));
        if (work_orders[i] == NULL)
            err(1, "impossible d'allouer un ordre de travail");
    }
    return work_orders;
}

int wait_for_work(int *work_order, int work_order_size)
{
    MPI_Status status;
    MPI_Recv(work_order, work_order_size, MPI_INTEGER, MPI_ANY_SOURCE,
            0, MPI_COMM_WORLD, &status);
    return status.MPI_SOURCE;
}

void send_work_to_proc(int worker, int *work_order, int work_order_size)
{
    MPI_Send(work_order, work_order_size, MPI_INTEGER, worker, 0, MPI_COMM_WORLD);
}

void launch_worker(const struct instance_t *instance, struct context_t *ctx,
        int work_order_size)
{
    int *work_order = malloc(work_order_size * sizeof(*work_order));
    if (work_order == NULL)
        err(1, "impossible d'allouer le buffer de distribution");

    work_order[0] = 0;          // distribution level
    work_order[1] = -1;         // last child at distrib_lvl (-1 means not defined)
    work_order[2] = rank - 1;   // first child at level 0
    do {
        int distrib_lvl = work_order[0];
        int last_child_at_distrib_lvl = work_order[1];

        // By default the first branch explored at each level is the first one
        memset(ctx->first_child, 0, instance->n_items * sizeof(*ctx->first_child));
        memcpy(ctx->first_child, work_order + 2, (distrib_lvl + 1) * sizeof(*work_order));

        print_work_order(work_order);
        solve(instance, ctx, work_order, work_order_size, nb_proc - 1,
                distrib_lvl, last_child_at_distrib_lvl);

        work_order[0] = -1;
        MPI_Send(work_order, work_order_size, MPI_INTEGER, 0, 0, MPI_COMM_WORLD);
        MPI_Recv(work_order, work_order_size, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } while(work_order[0] != -1);

}

void launch_master(int work_order_size)
{
    int nb_free_workers = 0;
    int *free_workers = malloc((nb_proc - 1) * sizeof(*free_workers));
    if (free_workers == NULL)
        err(1, "impossible d'allouer le tableau des travailleurs libres");

    int nb_work_orders = 0;
    int last_work_order = -1;
    int **work_orders = NULL;

    while (nb_free_workers < nb_proc - 1) {
        if (last_work_order + 1 >= nb_work_orders) {
            int new_nb_work_orders = nb_work_orders + INITIAL_NB_WORK_ORDERS;
            work_orders = allocate_work_orders(work_orders, nb_work_orders,
                    new_nb_work_orders, work_order_size);
            nb_work_orders = new_nb_work_orders;
        }

        int new_work_order = last_work_order + 1;
        /* Wait for work AND for done procs */
        int proc_src = wait_for_work(work_orders[new_work_order], work_order_size);
        if (work_orders[new_work_order][0] == -1) {
            free_workers[nb_free_workers++] = proc_src;
            DPRINTF("Master: proc [%d] is free\n", proc_src);
        } else {
            last_work_order = new_work_order;
            DPRINTF("Master: new work order from proc [%d] (total = %d)\n",
                    proc_src, last_work_order + 1);
        }

        /* Send remaining work to free workers */
        if (nb_free_workers > 0 && last_work_order >= 0) {
            send_work_to_proc(
                    free_workers[--nb_free_workers],
                    work_orders[last_work_order--],
                    work_order_size);
            DPRINTF("Master: sent work to proc [%d]\n", free_workers[nb_free_workers]);
        }
    }

    /* Announce to procs the work is done */
    work_orders[0][0] = -1;
    MPI_Request req;
    for (int proc = 1; proc < nb_proc; proc++)
        MPI_Isend(work_orders[0], work_order_size, MPI_INTEGER, proc, 0, MPI_COMM_WORLD, &req);
}

void launch_parallel(const struct instance_t *instance, struct context_t *ctx)
{
    /**
     * A work_order contains:
     * [distrib_level, num_children, child_num[0], ..., child_num[distrib_level]]
     * - num_children: number of children at the distrib_level
     */
    int work_order_size = 2 + instance->n_items;
    int root = 0;
    if (rank == root) {
        launch_master(work_order_size);
    } else {
        launch_worker(instance, ctx, work_order_size);
    }

    long long int total_solutions = 0;
    MPI_Reduce(&ctx->solutions, &total_solutions, 1, MPI_LONG_LONG_INT, MPI_SUM,
            root, MPI_COMM_WORLD);
    ctx->solutions = total_solutions;
}

struct context_t * backtracking_setup(const struct instance_t *instance)
{
    struct context_t *ctx = malloc(sizeof(*ctx));
    if (ctx == NULL)
        err(1, "impossible d'allouer un contexte");
    ctx->level = 0;
    ctx->nodes = 0;
    ctx->solutions = 0;
    ctx->next_work_share = 0;
    int n = instance->n_items;
    int m = instance->n_options;
    ctx->active_options = malloc(n * sizeof(*ctx->active_options));
    ctx->chosen_options = malloc(n * sizeof(*ctx->chosen_options));
    ctx->child_num = malloc(n * sizeof(*ctx->child_num));
    ctx->num_children = malloc(n * sizeof(*ctx->num_children));
    ctx->first_child = malloc(n * sizeof(*ctx->first_child));
    if (ctx->active_options == NULL || ctx->chosen_options == NULL
            || ctx->child_num == NULL || ctx->num_children == NULL
            || ctx->first_child == NULL)
        err(1, "impossible d'allouer le contexte");
    ctx->active_items = sparse_array_init(n);
    for (int item = 0; item < instance->n_primary; item++)
        sparse_array_add(ctx->active_items, item);

    for (int item = 0; item < n; item++)
        ctx->active_options[item] = sparse_array_init(m);
    for (int option = 0; option < m; option++)
        for (int k = instance->ptr[option]; k < instance->ptr[option + 1]; k++) {
            int item = instance->options[k];
            sparse_array_add(ctx->active_options[item], option);
        }

    return ctx;
}
