#include "datastructure.h"
#include <stdlib.h>
#include <err.h>


bool item_is_primary(const struct instance_t *instance, int item)
{
    return item < instance->n_primary;
}

struct sparse_array_t * sparse_array_init(int n)
{
    struct sparse_array_t *S = malloc(sizeof(*S));
    if (S == NULL)
        err(1, "impossible d'allouer un tableau creux");
    S->len = 0;
    S->capacity = n;
    S->p = malloc(n * sizeof(int));
    S->q = malloc(n * sizeof(int));
    if (S->p == NULL || S->q == NULL)
        err(1, "Impossible d'allouer p/q dans un tableau creux");

    // TODO: memset or parallel for loop
    for (int i = 0; i < n; i++)
        S->q[i] = n;           // initialement vide
    return S;
}

bool sparse_array_membership(const struct sparse_array_t *S, int x)
{
    return (S->q[x] < S->len);
}

bool sparse_array_empty(const struct sparse_array_t *S)
{
    return (S->len == 0);
}

void sparse_array_add(struct sparse_array_t *S, int x)
{
    int i = S->len;
    S->p[i] = x;
    S->q[x] = i;
    S->len = i + 1;
}

void sparse_array_remove(struct sparse_array_t *S, int x)
{
    int j = S->q[x];
    int n = S->len - 1;
    // échange p[j] et p[n]
    int y = S->p[n];
    S->p[n] = x;
    S->p[j] = y;
    // met q à jour
    S->q[x] = n;
    S->q[y] = j;
    S->len = n;
}

void sparse_array_unremove(struct sparse_array_t *S)
{
    S->len++;
}

void sparse_array_unadd(struct sparse_array_t *S)
{
    S->len--;
}



bool item_is_active(const struct context_t *ctx, int item)
{
    return sparse_array_membership(ctx->active_items, item);
}
