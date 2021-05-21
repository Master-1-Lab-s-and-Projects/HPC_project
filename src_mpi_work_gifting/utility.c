#include "utility.h"
#include <stdio.h>
#include <sys/time.h>
#include <err.h>


double wtime()
{
    struct timeval ts;
    gettimeofday(&ts, NULL);
    return (double) ts.tv_sec + ts.tv_usec / 1e6;
}

void print_option(const struct instance_t *instance, int option)
{
    if (instance->item_name == NULL)
        errx(1, "tentative d'affichage sans noms d'objet");
    for (int p = instance->ptr[option]; p < instance->ptr[option + 1]; p++) {
        int item = instance->options[p];
        printf("%s ", instance->item_name[item]);
    }
    printf("\n");
}

void print_sparse_array(struct sparse_array_t *arr)
{
    printf("[");
    for (int i=0; i < arr->len; i++)
        printf("%d%s", arr->p[i], (i == arr->len - 1)? "]\n": ", ");
}

void print_context(const struct context_t *ctx)
{
    printf("Context: \n");
    printf("* active_items: ");
    print_sparse_array(ctx->active_items);

    printf("* active_options: \n");
    for (int i=0; i < ctx->active_items->len; i++) {
        int item = ctx->active_items->p[i];
        printf("\t%d, ", item);
        print_sparse_array(ctx->active_options[item]);
    }
    printf("\n");
}

void print_work_order(int *work_order)
{
    int distrib_lvl = work_order[0];
    DPRINTF("Proc [%d]: New work order starts at level %d [", rank, distrib_lvl);
    for (int i = 0; i <= distrib_lvl; i++)
        DPRINTF("%d%s", work_order[i+2], (i == distrib_lvl) ? "" : ", ");
    DPRINTF(" - %d]\n", work_order[1]);
}

void progress_report(const struct context_t *ctx)
{
    double now = wtime();
    printf("Proc [%d] Exploré %lld noeuds, trouvé %lld solutions, temps écoulé %.1fs. ",
            rank, ctx->nodes, ctx->solutions, now - start);
    int i = 0;
    for (int k = 0; k < ctx->level; k++) {
        if (i > 44)
            break;
        int n = ctx->child_num[k];
        int m = ctx->num_children[k];
        if (m == 1)
            continue;
        printf("%c%c ", (n < 62) ? DIGITS[n] : '*', (m < 62) ? DIGITS[m] : '*');
        i++;
    }
    printf("\n"),
        next_report += report_delta;
}
