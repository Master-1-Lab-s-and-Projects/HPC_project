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
const long long WORK_SHARE_DELTA = 5e3; // partage son travail tous les ... noeuds
const int MAX_NB_OPTIONS_SHARED = 100;
int LOWEST_LEVEL_DISTRIBUTED = 0;
struct context_t **thread_ctxs = NULL;

double omp_time = 0.0;

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
    int d_lvl = 0; // highest level not completed
    if (LOWEST_LEVEL_DISTRIBUTED < lvl_limit)
        lvl_limit = LOWEST_LEVEL_DISTRIBUTED;

    if (d_lvl < lvl_limit && num_children[0] - child_num[0] <= nb_proc - 1)
        d_lvl++;

    while (d_lvl < lvl_limit && num_children[d_lvl] - child_num[d_lvl] <= 1)
        d_lvl++;

    if (num_children[d_lvl] - child_num[d_lvl] <= ((d_lvl == 0) ? nb_proc - 1 : 1))
        d_lvl = -1;

    return d_lvl;
}

// TODO
int get_first_child_to_share(const struct context_t *ctx, int distrib_lvl)
{
    int remaining_children = ctx->num_children[distrib_lvl] - ctx->child_num[distrib_lvl];
    assert(remaining_children > 1);
    if (distrib_lvl == 0) {
        DPRINTF("Proc [%d]-a: rem (%d)\n", rank, remaining_children);
        remaining_children -= remaining_children % (nb_proc - 1);
        DPRINTF("Proc [%d]-b: rem (%d)\n", rank, remaining_children);
        remaining_children /= nb_proc - 1;
        DPRINTF("Proc [%d]-c: rem (%d)\n", rank, remaining_children);
        remaining_children = ceil(((double) remaining_children)/ 2.0);
        DPRINTF("Proc [%d]-d: rem (%d)\n", rank, remaining_children);
        remaining_children *= nb_proc - 1;
        DPRINTF("Proc [%d]-e: rem (%d)\n", rank, remaining_children);
    } else {
        int children_to_share = min(remaining_children - 1, MAX_NB_OPTIONS_SHARED);
        remaining_children -= children_to_share;
    }

    return ctx->child_num[distrib_lvl] + remaining_children;
}

void share_work(struct context_t *ctx, int *work_order_buffer, int work_order_size)
{
    ctx->next_work_share += WORK_SHARE_DELTA;

    int distrib_lvl = get_distribution_lvl(ctx->child_num, ctx->num_children, ctx->level);
    if (distrib_lvl == -1) {
        DPRINTF("Proc [%d]: Did not share work\n", rank);
        return;
    }

    int increment = (distrib_lvl == 0) ? (nb_proc - 1) : 1;
    assert(ctx->child_num[distrib_lvl] + increment < ctx->num_children[distrib_lvl]);

    int first_child_shared = get_first_child_to_share(ctx, distrib_lvl);
    int last_child_shared = ctx->num_children[distrib_lvl];

    int n = ctx->active_items->capacity;
    int *shr_chosen_options = &work_order_buffer[3];
    int *shr_chosen_items   = &work_order_buffer[3 + n];
    int *shr_active_options = &work_order_buffer[3 + (2 * n)];

    memset(work_order_buffer, 0, (work_order_size) * sizeof(*work_order_buffer));
    work_order_buffer[0] = distrib_lvl;
    work_order_buffer[1] = first_child_shared;
    work_order_buffer[2] = last_child_shared;
    memcpy(shr_chosen_options, ctx->chosen_options, (distrib_lvl) * sizeof(*shr_chosen_options));
    memcpy(shr_chosen_items, ctx->chosen_items, (distrib_lvl + 1) * sizeof(*shr_chosen_items));

    DPRINTF("Proc [%02d]-1: SENT dlvl (%d)\n", rank, distrib_lvl);
    DPRINTF("Proc [%02d]-2: SENT chosen_options ", rank);
    print_array(ctx->chosen_options, 0, distrib_lvl);
    DPRINTF("Proc [%02d]-3: SENT chosen_items ", rank);
    print_array(ctx->chosen_items, 0, distrib_lvl + 1);

    if (distrib_lvl > 0) {
        int chosen_item_at_distrib_lvl = ctx->chosen_items[distrib_lvl];
        struct sparse_array_t *active_options = ctx->active_options[chosen_item_at_distrib_lvl];
        int *options_to_share = &active_options->p[first_child_shared];
        int nb_options_shared = last_child_shared - first_child_shared;
        memcpy(shr_active_options, options_to_share, nb_options_shared * sizeof(*shr_active_options));
        DPRINTF("Proc [%02d]-4: SENT shared_options ", rank);
        print_array(options_to_share, 0, nb_options_shared);
    }

    DPRINTF("Proc [%02d]-5: SENT ", rank);
    print_work_order(work_order_buffer);

    MPI_Send(work_order_buffer, work_order_size, MPI_INTEGER, 0, 0, MPI_COMM_WORLD);

    //DPRINTF("Proc [%d]: at level %d [%d], [%d, ..., %d[ --> [%d, ..., %d[\n",
    //        rank, distrib_lvl, (distrib_lvl > 0) ? ctx->chosen_options[distrib_lvl - 1] : -1,
    //        ctx->child_num[distrib_lvl], ctx->num_children[distrib_lvl],
    //        ctx->child_num[distrib_lvl], first_child_shared);
    ctx->num_children[distrib_lvl] = first_child_shared;
}

void solve2(const struct instance_t *instance, struct context_t *ctx)
{
    ctx->nodes++;
    if (ctx->nodes == next_report) {
        progress_report(ctx);
        printf("Proc [%d]: omp time %.3f\n", rank, omp_time);
    }
    if (sparse_array_empty(ctx->active_items)) {
        solution_found(instance, ctx);
        return;                         /* succès : plus d'objet actif */
    }
    int chosen_item = choose_next_item(ctx);
    struct sparse_array_t *active_options = ctx->active_options[chosen_item];
    if (sparse_array_empty(active_options))
        return;           /* échec : impossible de couvrir chosen_item */
    cover(instance, ctx, chosen_item);
    for (int k = 0; k < active_options->len; k++) {
        int option = active_options->p[k];
        choose_option(instance, ctx, option, chosen_item);
        solve2(instance, ctx);
        if (ctx->solutions >= max_solutions)
            return;
        unchoose_option(instance, ctx, option, chosen_item);
    }
    uncover(instance, ctx, chosen_item);                      /* backtrack */
}

void solve_omp(const struct instance_t *instance, struct context_t *ctx)
{
    ctx->nodes++;
    if (ctx->nodes == next_report) {
        progress_report(ctx);
        printf("Proc [%d]: omp time %.3f\n", rank, omp_time);
    }
    if (sparse_array_empty(ctx->active_items)) {
        solution_found(instance, ctx);
        return;                         /* succès : plus d'objet actif */
    }
    int chosen_item = choose_next_item(ctx);
    struct sparse_array_t *active_options = ctx->active_options[chosen_item];
    if (sparse_array_empty(active_options))
        return;           /* échec : impossible de couvrir chosen_item */
    cover(instance, ctx, chosen_item);

    long long solutions = ctx->solutions;
    long long nodes = ctx->nodes;
    #pragma omp parallel
    {
        int n_threads = omp_get_num_threads();

        struct context_t *pctx = thread_ctxs[omp_get_thread_num()];
        //printf("Proc [%d][%d]: pctx <%p> | ctx <%p>\n",
        //        rank, omp_get_thread_num(), pctx, ctx);
        copy_context(ctx, pctx);
        #pragma omp for schedule(runtime)
        for (int k = 0; k < active_options->len; k++) {
            int option = pctx->active_options[chosen_item]->p[k];
            choose_option(instance, pctx, option, chosen_item);
            solve2(instance, pctx);
            unchoose_option(instance, pctx, option, chosen_item);
        }

        #pragma omp single
        for (int i = 0; i < n_threads; i++) {
            ctx->solutions += thread_ctxs[i]->solutions - solutions;
            ctx->nodes += thread_ctxs[i]->nodes - nodes;
        }
    }
    uncover(instance, ctx, chosen_item);                      /* backtrack */
}

void solve(const struct instance_t *instance, struct context_t *ctx,
        int *work_order_buffer, int work_order_size, int loop_incr,
        int distrib_chosen_item, int first_child, int distrib_num_children)
{
    ctx->nodes++;
    if (ctx->nodes == next_report) {
        progress_report(ctx);
        printf("Proc [%d]: omp time %.3f\n", rank, omp_time);
    }
    if (sparse_array_empty(ctx->active_items)) {
        solution_found(instance, ctx);
        return;                         /* succès : plus d'objet actif */
    }
    int chosen_item = (distrib_chosen_item == -1)
        ? choose_next_item(ctx) : distrib_chosen_item;
    struct sparse_array_t *active_options = ctx->active_options[chosen_item];
    if (sparse_array_empty(active_options))
        return;           /* échec : impossible de couvrir chosen_item */
    cover(instance, ctx, chosen_item);

    ctx->num_children[ctx->level] = (distrib_num_children != -1)
        ? distrib_num_children : active_options->len;

    for (int k = first_child; k < ctx->num_children[ctx->level]; k += loop_incr) {
        int option = active_options->p[k];
        ctx->child_num[ctx->level] = k;
        if (ctx->nodes == ctx->next_work_share)
            share_work(ctx, work_order_buffer, work_order_size);
        choose_option(instance, ctx, option, chosen_item);
        if (ctx->level > LOWEST_LEVEL_DISTRIBUTED) {
            //printf("Proc [%d]: Entering OMP, lvl (%d)\n", rank, ctx->level);
            double start = wtime();
            solve_omp(instance, ctx);
            omp_time += wtime() - start;
            if (ctx->nodes >= ctx->next_work_share)
                ctx->next_work_share = ctx->nodes;
            if (ctx->nodes >= next_report)
                next_report = ctx->nodes;
        } else {
            solve(instance, ctx, work_order_buffer, work_order_size, 1, -1, 0, -1);
        }
        if (ctx->solutions >= max_solutions)
            return;
        unchoose_option(instance, ctx, option, chosen_item);
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
        const int *options_shared, int first_child_at_d_lvl, int last_child_at_d_lvl)
{
    int idx = 0;
    for (int i = first_child_at_d_lvl; i < last_child_at_d_lvl; i++) {
        int opt_a = active_options->p[i];
        int opt_b = options_shared[idx++];
        assert(sparse_array_membership(active_options, opt_b));
        int j = active_options->q[opt_b];

        active_options->p[i] = opt_b;
        active_options->q[opt_b] = i;

        active_options->p[j] = opt_a;
        active_options->q[opt_a] = j;
    }
}

void launch_worker(const struct instance_t *instance, struct context_t *ctx,
        int work_order_size)
{
    int *work_order = malloc(work_order_size * sizeof(*work_order));
    memset(work_order, 0, work_order_size * sizeof(*work_order));
    if (work_order == NULL)
        err(1, "impossible d'allouer le buffer de distribution");

    double idle_time = 0.0;
    double build_time = 0.0;
    int nb_work_orders_handled = 0;

    int n = instance->n_items;

    work_order[0] = 0;          // distribution level
    work_order[1] = rank - 1;   // first child at distrib_lvl
    work_order[2] = -1;         // last child at distrib_lvl (-1 means default)
    do {
        int distrib_lvl             = work_order[0];
        int first_child_at_d_lvl    = work_order[1];
        int last_child_at_d_lvl     = work_order[2];
        int *chosen_options         = &work_order[3];
        int *chosen_items           = &work_order[3 + n];
        int *options_shared         = &work_order[3 + (2 * n)];

        int loop_incr = (distrib_lvl == 0) ? (nb_proc - 1) : 1;
        int distrib_chosen_item = (distrib_lvl == 0) ? -1 : chosen_items[distrib_lvl];

        DPRINTF("Proc [%02d]: RCVD ....... - ", rank);
        print_work_order(work_order);

        double start_build = wtime();
        build_context(instance, ctx, chosen_items, chosen_options, distrib_lvl);
        if (ctx->level > 0)
            reorder_active_options(ctx->active_options[distrib_chosen_item],
                    options_shared, first_child_at_d_lvl, last_child_at_d_lvl);
        build_time += wtime() - start_build;

        solve(instance, ctx, work_order, work_order_size, loop_incr,
                distrib_chosen_item, first_child_at_d_lvl, last_child_at_d_lvl);
        nb_work_orders_handled++;

        double start_reset = wtime();
        reset_context(instance, ctx, ctx->chosen_items, ctx->chosen_options, distrib_lvl);
        build_time += wtime() - start_reset;

        double start_idle = wtime();

        work_order[0] = -1;
        MPI_Send(work_order, work_order_size, MPI_INTEGER, 0, 0, MPI_COMM_WORLD);
        MPI_Recv(work_order, work_order_size, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // TODO: count waiting time when no more work
        if (work_order[0] != -1)
            idle_time += wtime() - start_idle;

    } while(work_order[0] != -1);
            DPRINTF("Proc [%d]: solutions (%lld), nb_work_orders_handled (%d)\n", rank, ctx->solutions, nb_work_orders_handled);
            DPRINTF("Proc [%d]: idle_time = %.3fs, build_time = %.3fs\n", rank, idle_time, build_time);
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

    int max_work_orders = 0;
    int deepest_level_shrd = 0;
    while (nb_free_workers < nb_proc - 1) {
        if (nb_work_orders >= work_orders_capacity)
            work_orders = allocate_work_orders(work_orders, &work_orders_capacity, work_order_size);

        /* Wait for work AND for done procs */
        double start_idle = wtime();
        int proc_src = wait_for_work(work_orders[nb_work_orders], work_order_size);
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
                    work_orders[--nb_work_orders],
                    work_order_size);
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
        MPI_Isend(work_orders[0], work_order_size, MPI_INTEGER, proc, 0, MPI_COMM_WORLD, &req);
    DPRINTF("Master: idle_time = %.3fs, Max work orders at the same time = %d, Deepest lvl shrd = %d\n",
            idle_time, max_work_orders, deepest_level_shrd);

    free(free_workers);
    for (int i = 0; i < work_orders_capacity; i++)
        free(work_orders[i]);
    free(work_orders);
}

void setup_omp(const struct context_t *ctx)
{
    int n_max_threads = omp_get_max_threads();
    thread_ctxs = malloc(n_max_threads * sizeof(*thread_ctxs));
    if (thread_ctxs == NULL)
        err(1, "Impossible d'allouer le tableau de context pour les threads");

    #pragma omp parallel for
    for (int i = 0; i < n_max_threads; i++)
        thread_ctxs[i] = new_copy_context(ctx);
}

void launch_parallel(const struct instance_t *instance, struct context_t *ctx)
{
    /**
     * A work_order wo contains:
     * - wo[0]: distrib_lvl
     * - wo[1]: first_child_shared
     * - wo[2]: last_child_shared
     * - wo[3]: chosen_options (up to [n] options)
     * - wo[3 + n]: chosen_items (up to [n] items)
     * - wo[3 + (2 * n)]: options_shared (up to [MAX_NB_OPTIONS_SHARED] options)
     */
    int n = instance->n_items;
    int work_order_size = 3 + (2 * n) + MAX_NB_OPTIONS_SHARED;
    LOWEST_LEVEL_DISTRIBUTED = ceil(n * 0.1); // TODO: Find a name for this constant/parameter

    setup_omp(ctx);

    if (nb_proc < 3)
        ctx->next_work_share = -1;

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
    ctx->chosen_items   = malloc(n * sizeof(*ctx->chosen_options));
    ctx->child_num      = malloc(n * sizeof(*ctx->child_num));
    ctx->num_children   = malloc(n * sizeof(*ctx->num_children));
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
