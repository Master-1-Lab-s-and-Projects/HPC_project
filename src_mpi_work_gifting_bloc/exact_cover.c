#include "exact_cover.h"
#include "utility.h"
#include <stdio.h>
#include <stdlib.h>
#include <err.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include <mpi.h>
#include <omp.h>

#define INCR_WORK_ORDERS_CAPACITY 2 * nb_proc
const long long WORK_SHARE_DELTA = 5e5; // partage son travail tous les ... noeuds
const int MAX_NB_OPTIONS_SHARED = 500;
int MAX_LEVEL_SHARED = -1;
int WORK_ORDER_SIZE = -1;

int first_level = -1;

struct context_t **thr_ctxs = NULL;     // context used by each thread
int **thr_work_order_buffers = NULL;    //  work order buffers for each thread
struct sparse_array_t *fst_lvl_active_options = NULL;
omp_lock_t fst_lvl_lock;


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
    ctx->chosen_items[ctx->level] = item;
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
    int d_lvl = first_level + 1; // highest level not completed
    if (MAX_LEVEL_SHARED < lvl_limit)
        lvl_limit = MAX_LEVEL_SHARED;

    while (d_lvl < lvl_limit && num_children[d_lvl] - child_num[d_lvl] <= 1)
        d_lvl++;

    if (d_lvl > lvl_limit || num_children[d_lvl] - child_num[d_lvl] <= 1)
        d_lvl = -1;

    return d_lvl;
}

// TODO
int get_first_child_to_share(const struct context_t *ctx, int distrib_lvl)
{
    if (distrib_lvl == first_level)
        return 0;
    int remaining_children = ctx->num_children[distrib_lvl] - ctx->child_num[distrib_lvl];
    int half_remain = ceil(((double) remaining_children)/2.0);
    int children_to_share = min(half_remain, MAX_NB_OPTIONS_SHARED);
    remaining_children -= children_to_share;

    return ctx->child_num[distrib_lvl] + remaining_children;
}

void share_work(const struct instance_t *instance, struct context_t *ctx,
        int *work_order_buffer)
{
    ctx->next_work_share += WORK_SHARE_DELTA;

    omp_set_lock(&fst_lvl_lock);

    int distrib_lvl = first_level;
    if (sparse_array_empty(fst_lvl_active_options))
        distrib_lvl = get_distribution_lvl(ctx->child_num, ctx->num_children, ctx->level);
    if (distrib_lvl == -1) {
        DPRINTF("Proc [%d]: Did not share work\n", rank);
        omp_unset_lock(&fst_lvl_lock);
        return;
    }

    int first_child_shared = get_first_child_to_share(ctx, distrib_lvl);
    int last_child_shared = (distrib_lvl == first_level)
        ? min(fst_lvl_active_options->len, MAX_NB_OPTIONS_SHARED)
        :ctx->num_children[distrib_lvl];
    int nb_options_shared = last_child_shared - first_child_shared;

    int n = instance->n_items;
    int *shr_chosen_options = &work_order_buffer[2];
    int *shr_chosen_items   = &work_order_buffer[2 + n];
    int *shr_active_options = &work_order_buffer[2 + (2 * n)];

    memset(work_order_buffer, 0, (WORK_ORDER_SIZE) * sizeof(*work_order_buffer));
    work_order_buffer[0] = distrib_lvl;
    work_order_buffer[1] = nb_options_shared;
    memcpy(shr_chosen_options, ctx->chosen_options, (distrib_lvl) * sizeof(*shr_chosen_options));
    memcpy(shr_chosen_items, ctx->chosen_items, (distrib_lvl + 1) * sizeof(*shr_chosen_items));

    int *options_to_share = NULL;
    if (distrib_lvl == first_level) {
        options_to_share = fst_lvl_active_options->p;
    } else {
        int chosen_item = ctx->chosen_items[distrib_lvl];
        struct sparse_array_t *active_opts = ctx->active_options[chosen_item];
        options_to_share = &active_opts->p[first_child_shared];
    }

    memcpy(shr_active_options, options_to_share, nb_options_shared * sizeof(*shr_active_options));
    if (distrib_lvl == first_level) {
        int len_before = fst_lvl_active_options->len;
        for (int i = 0; i < nb_options_shared; i++) {
            assert(sparse_array_membership(fst_lvl_active_options, shr_active_options[i]));
            sparse_array_remove(fst_lvl_active_options, shr_active_options[i]);
        }
        assert(len_before - nb_options_shared == fst_lvl_active_options->len);
    }
    omp_unset_lock(&fst_lvl_lock);

    MPI_Send(work_order_buffer, WORK_ORDER_SIZE, MPI_INTEGER, 0, 0, MPI_COMM_WORLD);
    ctx->num_children[distrib_lvl] = first_child_shared;
}

void solve2(const struct instance_t *instance, struct context_t *ctx,
        int *work_order_buffer)
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
    for (int k = 0; k < ctx->num_children[ctx->level]; k++) {
        int option = active_options->p[k];
        ctx->child_num[ctx->level] = k;
        if (ctx->nodes == ctx->next_work_share)
            share_work(instance, ctx, work_order_buffer);
        choose_option(instance, ctx, option, chosen_item);
        solve2(instance, ctx, work_order_buffer);
        if (ctx->solutions >= max_solutions)
            return;
        unchoose_option(instance, ctx, option, chosen_item);
    }
    uncover(instance, ctx, chosen_item);                      /* backtrack */
}

void solve(const struct instance_t *instance, struct context_t *ctx,
        int chosen_item, int num_children)
{
    const struct sparse_array_t *active_options = ctx->active_options[chosen_item];
    if (sparse_array_empty(active_options))
        return;           /* échec : impossible de couvrir chosen_item */
    cover(instance, ctx, chosen_item);

    long long nodes_before = ctx->nodes;
    long long solutions_before = ctx->solutions;
    copy_sparse_array(active_options, fst_lvl_active_options);
    fst_lvl_active_options->len = num_children;
    assert(first_level == ctx->level);
    #pragma omp parallel
    {
        int thr_num = omp_get_thread_num();
        int *work_order_buffer = thr_work_order_buffers[thr_num];

        struct context_t *pctx = thr_ctxs[thr_num];
        copy_context(ctx, pctx);
        pctx->num_children[pctx->level] = -1;

        #pragma omp for schedule(dynamic, 1)
        for (int k = 0; k < num_children; k++) {
            int option = active_options->p[k];
            bool was_shared = true;

            omp_set_lock(&fst_lvl_lock);
            if (sparse_array_membership(fst_lvl_active_options, option)) {
                sparse_array_remove(fst_lvl_active_options, option);
                was_shared = false;
            }
            omp_unset_lock(&fst_lvl_lock);

            if (!was_shared) {
                if (pctx->nodes == pctx->next_work_share)
                    share_work(instance, pctx, work_order_buffer);
                choose_option(instance, pctx, option, chosen_item);
                solve2(instance, pctx, work_order_buffer);
                unchoose_option(instance, pctx, option, chosen_item);
            }
        }
        long long nodes_seen = pctx->nodes - nodes_before;
        long long solutions_seen = pctx->solutions - solutions_before;

        #pragma omp atomic
        ctx->nodes += nodes_seen;
        #pragma omp atomic
        ctx->solutions += solutions_seen;
    }
    assert(fst_lvl_active_options->len == 0);

    ctx->next_work_share = ctx->nodes - (ctx->nodes % WORK_SHARE_DELTA);
    ctx->next_work_share += WORK_SHARE_DELTA;
    uncover(instance, ctx, chosen_item);                      /* backtrack */
}

int ** allocate_work_orders(int **work_orders, int *work_orders_capacity)
{
    int old_capacity = *work_orders_capacity;
    int new_capacity = old_capacity + INCR_WORK_ORDERS_CAPACITY;
    work_orders = realloc(work_orders, new_capacity * sizeof(*work_orders));
    if (work_orders == NULL)
        err(1, "impossible d'allouer le tableau des ordres de travail");

    for (int i = old_capacity; i < new_capacity; i++) {
        work_orders[i] = malloc(WORK_ORDER_SIZE * sizeof(*work_orders[i]));
        if (work_orders[i] == NULL)
            err(1, "impossible d'allouer un ordre de travail");
    }

    *work_orders_capacity = new_capacity;
    return work_orders;
}

int wait_for_work(int *work_order)
{
    MPI_Status status;
    MPI_Recv(work_order, WORK_ORDER_SIZE, MPI_INTEGER, MPI_ANY_SOURCE,
            0, MPI_COMM_WORLD, &status);
    return status.MPI_SOURCE;
}

void send_work_to_proc(int worker, int *work_order)
{
    MPI_Send(work_order, WORK_ORDER_SIZE, MPI_INTEGER, worker, 0, MPI_COMM_WORLD);
}

void build_context(const struct instance_t *instance, struct context_t *ctx,
        const int *chosen_items, const int *chosen_options, int distrib_lvl)
{
    for (int i = 0; i < distrib_lvl; i++) {
        ctx->child_num[i] = 0;
        ctx->num_children[i] = 1;
        cover(instance, ctx, chosen_items[i]);
        choose_option(instance, ctx, chosen_options[i], chosen_items[i]);
    }
}

void reset_context(const struct instance_t *instance, struct context_t *ctx,
        const int *chosen_items, const int *chosen_options, int distrib_lvl)
{
    for (int i = distrib_lvl - 1; i >= 0; i--) {
        unchoose_option(instance, ctx, chosen_options[i], chosen_items[i]);
        uncover(instance, ctx, chosen_items[i]);
    }
}

void reorder_active_options(struct sparse_array_t *active_options,
        const int *options_shared, int nb_options)
{
    for (int i = 0; i < nb_options; i++) {
        int opt_a = active_options->p[i];
        int opt_b = options_shared[i];
        assert(sparse_array_membership(active_options, opt_b));
        int j = active_options->q[opt_b];

        active_options->p[i] = opt_b;
        active_options->q[opt_b] = i;

        active_options->p[j] = opt_a;
        active_options->q[opt_a] = j;
    }
}

void launch_worker(const struct instance_t *instance, struct context_t *ctx)
{
    int *work_order = malloc(WORK_ORDER_SIZE * sizeof(*work_order));
    memset(work_order, 0, WORK_ORDER_SIZE * sizeof(*work_order));
    if (work_order == NULL)
        err(1, "impossible d'allouer le buffer de distribution");

    double idle_time = 0.0;
    double build_time = 0.0;
    int nb_work_orders_handled = 0;

    int n = instance->n_items;

    int chosen_item = choose_next_item(ctx);
    struct sparse_array_t *active_options = ctx->active_options[chosen_item];
    int nb_loops = active_options->len;
    int chunk = nb_loops / (nb_proc - 1);
    ctx->nodes++;
    if (sparse_array_empty(ctx->active_items)) {
        solution_found(instance, ctx);
        return;                         /* succès : plus d'objet actif */
    }
    int fst_opt_idx = chunk * (rank - 1);
    int lst_opt_idx = (rank == nb_proc - 1) ? nb_loops : fst_opt_idx + chunk;
    int nb_opts = lst_opt_idx - fst_opt_idx;

    copy_sparse_array(active_options, fst_lvl_active_options);
    reorder_active_options(active_options, fst_lvl_active_options->p + fst_opt_idx, nb_opts);

    work_order[0] = 0;          // distribution level
    work_order[1] = nb_opts;
    work_order[2 + n] = chosen_item;
    bool first_work = true;
    do {
        first_level = work_order[0];
        int distrib_lvl     = work_order[0];
        int nb_options      = work_order[1];
        int *chosen_options = &work_order[2];
        int *chosen_items   = &work_order[2 + n];
        int *options_shared = &work_order[2 + (2 * n)];

        int d_chosen_item = chosen_items[distrib_lvl];

        DPRINTF("Proc [%d]: RCVD next (%lld), curr (%lld)\n", rank, ctx->next_work_share, ctx->nodes);
        print_work_order(work_order);

        double start_build = wtime();
        build_context(instance, ctx, chosen_items, chosen_options, distrib_lvl);
        if (!first_work)
            reorder_active_options(ctx->active_options[d_chosen_item], options_shared, nb_options);
        build_time += wtime() - start_build;

        solve(instance, ctx, d_chosen_item, nb_options);
        nb_work_orders_handled++;
        first_work = false;

        double start_reset = wtime();
        reset_context(instance, ctx, ctx->chosen_items, ctx->chosen_options, distrib_lvl);
        build_time += wtime() - start_reset;

        double start_idle = wtime();

        work_order[0] = -1;
        MPI_Send(work_order, WORK_ORDER_SIZE, MPI_INTEGER, 0, 0, MPI_COMM_WORLD);
        MPI_Recv(work_order, WORK_ORDER_SIZE, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // TODO: count waiting time when no more work
        if (work_order[0] != -1)
            idle_time += wtime() - start_idle;

    } while(work_order[0] != -1);
    printf("Proc [%d]: solutions (%lld), nb_work_orders_handled (%d)\n", rank, ctx->solutions, nb_work_orders_handled);
    printf("Proc [%d]: idle_time = %.3fs, build_time = %.3fs\n", rank, idle_time, build_time);
    free(work_order);
}

void swap_work_orders(int **work_order_a, int **work_order_b)
{
    int *tmp = *work_order_a;
    *work_order_a = *work_order_b;
    *work_order_b = tmp;
}

void swap_free_workers(int *worker_a, int *worker_b)
{
    int tmp = *worker_a;
    *worker_a = *worker_b;
    *worker_b = tmp;
}

void launch_master()
{
    int nb_free_workers = 0;
    int *free_workers = malloc((nb_proc - 1) * sizeof(*free_workers));
    if (free_workers == NULL)
        err(1, "impossible d'allouer le tableau des travailleurs libres");

    int work_orders_capacity = 0;
    int nb_work_orders = 0;
    int **work_orders = NULL;

    double idle_time = 0.0;

    int max_work_orders = 0;
    int deepest_level_shrd = 0;
    while (nb_free_workers < nb_proc - 1) {
        if (nb_work_orders >= work_orders_capacity)
            work_orders = allocate_work_orders(work_orders, &work_orders_capacity);

        /* Wait for work AND for done procs */
        double start_idle = wtime();
        int proc_src = wait_for_work(work_orders[nb_work_orders]);
        idle_time += wtime() - start_idle;

        if (work_orders[nb_work_orders][0] == -1) {
            DPRINTF("Master: RCVD free proc n°%d\n", proc_src);
            free_workers[nb_free_workers++] = proc_src;
            //if (nb_free_workers > 0)
            //    swap_free_workers(&free_workers[0], &free_workers[nb_free_workers - 1]);
        } else {
            nb_work_orders++;
            DPRINTF("Master: RCVD work order from proc n°%d: ", proc_src);
            print_work_order(work_orders[nb_work_orders - 1]);
            //if (nb_work_orders > 1)
            //    swap_work_orders(&work_orders[0], &work_orders[nb_work_orders - 1]);
            if (work_orders[nb_work_orders - 1][0] > deepest_level_shrd)
                deepest_level_shrd = work_orders[nb_work_orders - 1][0];
            if (nb_work_orders > max_work_orders)
                max_work_orders = nb_work_orders;
        }

        /* Send remaining work to free workers */
        if (nb_free_workers > 0 && nb_work_orders > 0) {
            DPRINTF("Master: SENT to proc %02d - ", free_workers[nb_free_workers - 1]);
            print_work_order(work_orders[nb_work_orders - 1]);
            send_work_to_proc(
                    free_workers[--nb_free_workers],
                    work_orders[--nb_work_orders]);
            assert(nb_free_workers >= 0);
            assert(nb_work_orders >= 0);
        }

        DPRINTF("Master: work orders (%d) vvv, free procs (%d): ",
                nb_work_orders, nb_free_workers);
        print_array(free_workers, 0, nb_free_workers);
        for (int i = 0; i < nb_work_orders; i++) {
            DPRINTF("  * ");
            print_work_order(work_orders[i]);
        }
    }
    assert(nb_free_workers == nb_proc - 1);
    assert(nb_work_orders == 0);

    /* Announce to procs the work is done */
    work_orders[0][0] = -1;
    MPI_Request req;
    for (int proc = 1; proc < nb_proc; proc++)
        MPI_Isend(work_orders[0], WORK_ORDER_SIZE, MPI_INTEGER, proc, 0, MPI_COMM_WORLD, &req);
    printf("Master: idle_time = %.3fs, Max work orders at the same time = %d, Deepest lvl shrd = %d\n",
            idle_time, max_work_orders, deepest_level_shrd);

    free(free_workers);
    for (int i = 0; i < work_orders_capacity; i++)
        free(work_orders[i]);
    free(work_orders);
}

void setup_omp(const struct context_t *ctx, int n_options)
{
    int n_max_threads = omp_get_max_threads();
    thr_ctxs = malloc(n_max_threads * sizeof(*thr_ctxs));
    if (thr_ctxs == NULL)
        err(1, "Impossible d'allouer le tableau de contexts pour les threads");

    //#pragma omp parallel for
    for (int i = 0; i < n_max_threads; i++)
        thr_ctxs[i] = new_copy_context(ctx);

    thr_work_order_buffers = malloc(n_max_threads * sizeof(*thr_work_order_buffers));
    if (thr_work_order_buffers == NULL)
        err(1, "Impossible d'allouer le tableau de work orders pour les threads");

    //#pragma omp parallel for
    for (int i = 0; i < n_max_threads; i++) {
        thr_work_order_buffers[i] = malloc(WORK_ORDER_SIZE * sizeof(*thr_work_order_buffers[i]));
        if (thr_work_order_buffers[i] == NULL)
            err(1, "Impossible d'allouer le work order pour le thread");
    }

    fst_lvl_active_options = sparse_array_init(n_options);
    omp_init_lock(&fst_lvl_lock);
}

void launch_parallel(const struct instance_t *instance, struct context_t *ctx)
{
    /**
     * A work_order wo contains:
     * - wo[0]: distrib_lvl
     * - wo[1]: last_child_shared
     * - wo[2]: chosen_options (up to [n] options)
     * - wo[2 + n]: chosen_items (up to [n] items)
     * - wo[2 + (2 * n)]: options_shared (up to [MAX_NB_OPTIONS_SHARED] options)
     */
    int n = instance->n_items;
    WORK_ORDER_SIZE = 2 + (2 * n) + MAX_NB_OPTIONS_SHARED;
    MAX_LEVEL_SHARED = ceil(n * 0.4); // TODO: Find a name for this constant/parameter

    setup_omp(ctx, instance->n_options);

    int root = 0;
    if (rank == root)
        launch_master();
    else
        launch_worker(instance, ctx);

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
    ctx->chosen_items = malloc(n * sizeof(*ctx->chosen_options));
    ctx->child_num = malloc(n * sizeof(*ctx->child_num));
    ctx->num_children = malloc(n * sizeof(*ctx->num_children));
    if (ctx->active_options == NULL || ctx->chosen_options == NULL
            || ctx->child_num == NULL || ctx->num_children == NULL)
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
