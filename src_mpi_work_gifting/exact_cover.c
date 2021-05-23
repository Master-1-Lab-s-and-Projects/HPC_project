#include "exact_cover.h"
#include "utility.h"
#include <stdio.h>
#include <stdlib.h>
#include <err.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include <mpi.h>

#define INCR_WORK_ORDERS_CAPACITY 2 * nb_proc
const long long WORK_SHARE_DELTA = 1e5; // partage son travail tous les ... noeuds
const int MAX_LEVEL_SHARED = 1;


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
        int option, int chosen_item, bool was_chosen)
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
        int option, int chosen_item, bool was_chosen)
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
int get_distribution_lvl(const int *child_num, const int* num_children, int lvl_limit)
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
    assert(ctx->child_num[distrib_lvl] < first_child_shared
            && first_child_shared < last_child_shared);

    memcpy(&work_order_buffer[3], ctx->chosen_options, (distrib_lvl) * sizeof(*ctx->chosen_options));
    work_order_buffer[0] = distrib_lvl;
    work_order_buffer[1] = first_child_shared;
    work_order_buffer[2] = last_child_shared;
    MPI_Send(work_order_buffer, work_order_size, MPI_INTEGER, 0, 0, MPI_COMM_WORLD);


    DPRINTF("Proc [%d]: at level %d [%d], [%d, ..., %d[ --> [%d, ..., %d[\n",
            rank, distrib_lvl, (distrib_lvl == 1) ? ctx->chosen_options[0] : -1,
            ctx->child_num[distrib_lvl], ctx->num_children[distrib_lvl],
            ctx->child_num[distrib_lvl], first_child_shared);
    ctx->num_children[distrib_lvl] = first_child_shared;

}

void solve(const struct instance_t *instance, struct context_t *ctx,
        int *work_order_buffer, int work_order_size, int loop_incr,
        int distrib_lvl, int distrib_num_children)
{
    assert(ctx->level >= distrib_lvl);
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

    if (ctx->level == -1) {
        DPRINTF("Proc [%d]: co (%d) ci (%d), ", rank, ctx->chosen_options[0], chosen_item);
        print_sparse_array(ctx->active_items);
    }
    cover(instance, ctx, chosen_item);

    ctx->num_children[ctx->level] = (distrib_num_children != -1)
        ? distrib_num_children : active_options->len;

    int old_opt = 0;
    if (ctx->level == 1)
        old_opt = ctx->chosen_options[0];
    for (int k = ctx->first_child[ctx->level]; k < ctx->num_children[ctx->level]; k += loop_incr) {
        int option = active_options->p[k];
        ctx->child_num[ctx->level] = k;
        choose_option(instance, ctx, option, chosen_item, true);
        solve(instance, ctx, work_order_buffer, work_order_size, 1, distrib_lvl, -1);

        if (ctx->solutions >= max_solutions)
            return;
        unchoose_option(instance, ctx, option, chosen_item, true);
    }
    if (ctx->level == -1) {
        DPRINTF("%d, %d ... %d (Proc [%d])\n",
                ctx->chosen_options[0],
                ctx->first_child[1], ctx->num_children[1], rank);
        assert(old_opt == ctx->chosen_options[0]);
    }

    uncover(instance, ctx, chosen_item);                      /* backtrack */
}

int ** allocate_work_orders(int **work_orders, int *work_orders_capacity,
        int work_order_size)
{
    int old_capacity = *work_orders_capacity;
    int new_capacity = old_capacity + INCR_WORK_ORDERS_CAPACITY;
    work_orders = realloc(work_orders, new_capacity * sizeof(*work_orders));
    if (work_orders == NULL)
        err(1, "impossible d'allouer le tableau des ordres de travail");

    for (int i = old_capacity; i < new_capacity; i++) {
        work_orders[i] = malloc(work_order_size * sizeof(*work_orders[i]));
        if (work_orders[i] == NULL)
            err(1, "impossible d'allouer un ordre de travail");
    }

    *work_orders_capacity = new_capacity;
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

void build_context(const struct instance_t *instance,
        struct context_t *ctx, const int *chosen_options, int n_options)
{
    for (int i = 0; i < n_options; i++)
        choose_option(instance, ctx, chosen_options[i], -1, false);
}

void reset_context(const struct instance_t *instance,
        struct context_t *ctx, const int *chosen_options, int n_options)
{
    for (int i = n_options - 1; i >= 0; i--)
        unchoose_option(instance, ctx, chosen_options[i], -1, false);
}

void new_sparse_array(const struct sparse_array_t *src, struct sparse_array_t *dst)
{
    dst->capacity = src->capacity;
    dst->p = malloc(dst->capacity * sizeof(*dst->p));
    dst->q = malloc(dst->capacity * sizeof(*dst->q));
}

void copy_sparse_array(const struct sparse_array_t *src, struct sparse_array_t *dst)
{
    assert(dst->capacity == src->capacity);
    dst->len = src->len;
    memcpy(dst->p, src->p, dst->capacity * sizeof(*dst->p));
    memcpy(dst->q, src->q, dst->capacity * sizeof(*dst->q));
}

void launch_worker(const struct instance_t *instance, struct context_t *ctx,
        int work_order_size)
{
    int *work_order = malloc(work_order_size * sizeof(*work_order));
    if (work_order == NULL)
        err(1, "impossible d'allouer le buffer de distribution");

    double idle_time = 0.0;

    int n = instance->n_items;
    struct sparse_array_t active_items_sav = {0};
    new_sparse_array(ctx->active_items, &active_items_sav);
    copy_sparse_array(ctx->active_items, &active_items_sav);

    struct sparse_array_t *active_options_sav = malloc(n * sizeof(*ctx->active_options));
    for (int i = 0; i < n; i++) {
        new_sparse_array(ctx->active_options[0], &active_options_sav[0]);
        copy_sparse_array(ctx->active_options[0], &active_options_sav[0]);
    }



    work_order[0] = 0;          // distribution level
    work_order[1] = rank - 1;   // first child at distrib_lvl
    work_order[2] = -1;         // last child at distrib_lvl (-1 means default)
    do {
        int distrib_lvl = work_order[0];
        int first_child_at_d_lvl = work_order[1];
        int last_child_at_d_lvl = work_order[2];

        // By default the first branch explored at each level is the first one
        memset(ctx->first_child, 0, instance->n_items * sizeof(*ctx->first_child));
        ctx->first_child[distrib_lvl] = first_child_at_d_lvl;

        memcpy(ctx->chosen_options, &work_order[3], distrib_lvl * sizeof(ctx->chosen_options));

        int loop_incr = (distrib_lvl == 0) ? (nb_proc - 1) : 1;

        build_context(instance, ctx, ctx->chosen_options, distrib_lvl);
        solve(instance, ctx, work_order, work_order_size, loop_incr,
                distrib_lvl, last_child_at_d_lvl);
        reset_context(instance, ctx, ctx->chosen_options, distrib_lvl);

        copy_sparse_array(&active_items_sav, ctx->active_items);
        for (int i = 0; i < n; i++)
            copy_sparse_array(&active_options_sav[0], ctx->active_options[0]);
        ctx->level = 0;

        work_order[0] = -1;
        double start_idle = wtime();
        MPI_Send(work_order, work_order_size, MPI_INTEGER, 0, 0, MPI_COMM_WORLD);
        MPI_Recv(work_order, work_order_size, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (work_order[0] != -1)
            idle_time += wtime() - start_idle;
    } while(work_order[0] != -1);

    printf("Proc [%d]: idle_time = %.3fs\n", rank, idle_time);
    free(work_order);
}

void launch_master(int work_order_size)
{
    int nb_free_workers = 0;
    int *free_workers = malloc((nb_proc - 1) * sizeof(*free_workers));
    if (free_workers == NULL)
        err(1, "impossible d'allouer le tableau des travailleurs libres");

    int work_orders_capacity = 0;
    int nb_work_orders = 0;
    int **work_orders = NULL;

    double idle_time = 0.0;

    while (nb_free_workers < nb_proc - 1) {
        if (nb_work_orders >= work_orders_capacity)
            work_orders = allocate_work_orders(work_orders, &work_orders_capacity, work_order_size);

        double start_idle = wtime();
        //int new_work_order = nb_work_orders;
        /* Wait for work AND for done procs */
        int proc_src = wait_for_work(work_orders[nb_work_orders], work_order_size);
        idle_time += wtime() - start_idle;
        if (work_orders[nb_work_orders][0] == -1) {
            free_workers[nb_free_workers++] = proc_src;
            //DPRINTF("Master: proc [%d] is free\n", proc_src);
        } else {
            //if (nb_work_orders > 0) {
            //    int *last_work_order = work_orders[0];
            //    work_orders[0] = work_orders[nb_work_orders];
            //    work_orders[nb_work_orders] = last_work_order;
            //}
            nb_work_orders++;
            //DPRINTF("Master: new work order from proc [%d] (total = %d): ",
            //        proc_src, nb_work_orders);
            //print_work_order(work_orders[nb_work_orders - 1]);
        }

        /* Send remaining work to free workers */
        if (nb_free_workers > 0 && nb_work_orders > 0) {
            send_work_to_proc(
                    free_workers[--nb_free_workers],
                    work_orders[--nb_work_orders],
                    work_order_size);
            //DPRINTF("Master: sent work to proc [%d] (total = %d)\n",
            //        free_workers[nb_free_workers], nb_work_orders);
        }
        //for (int i = 0; i < nb_work_orders; i++) {
        //    DPRINTF("Master: (%d) ", i);
        //    print_work_order(work_orders[i]);
        //}
    }

    /* Announce to procs the work is done */
    work_orders[0][0] = -1;
    MPI_Request req;
    for (int proc = 1; proc < nb_proc; proc++)
        MPI_Isend(work_orders[0], work_order_size, MPI_INTEGER, proc, 0, MPI_COMM_WORLD, &req);
    printf("Master: idle_time = %.3fs\n", idle_time);

    for (int i = 0; i < work_orders_capacity; i++)
        free(work_orders[i]);
    free(work_orders);
    free(free_workers);
}

void launch_parallel(const struct instance_t *instance, struct context_t *ctx)
{
    /**
     * A work_order contains:
     * [distrib_lvl, first_child, last_child, option[0], ..., child_num[distrib_lvl]]
     * - first_child: first child to visit at distrib_lvl
     * - last_child: last child (excluded) to visit at distrib_lvl
     */
    int work_order_size = 3 + instance->n_items;
    int root = 0;
    if (rank == root)
        launch_master(work_order_size);
    else
        launch_worker(instance, ctx, work_order_size);

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
    ctx->next_work_share = WORK_SHARE_DELTA;
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
