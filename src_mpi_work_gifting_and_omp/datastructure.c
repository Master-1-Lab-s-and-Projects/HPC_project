#include "datastructure.h"
#include <stdlib.h>
#include <err.h>
#include <assert.h>
#include <string.h>

int * new_copy_int_array(const int *arr, int len)
{
    int *cpy = malloc(len * sizeof(*cpy));
    if (cpy == NULL)
        err(1, "impossible d'allouer un tableau d'entiers");
    memcpy(cpy, arr, len * sizeof(*cpy));
    return cpy;
}

void copy_int_array(const int *src, int *dst, int len)
{
    memcpy(dst, src, len * sizeof(*src));
}

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


struct sparse_array_t * new_copy_sparse_array(const struct sparse_array_t *S)
{
    struct sparse_array_t *cpy = malloc(sizeof(*cpy));
    if (cpy == NULL)
        err(1, "impossible d'allouer un tableau creux");
    cpy->capacity = S->capacity;
    cpy->len = S->len;
    cpy->p = new_copy_int_array(S->p, S->capacity);
    cpy->q = new_copy_int_array(S->q, S->capacity);
    return cpy;
}

void copy_sparse_array(const struct sparse_array_t *src,
        struct sparse_array_t *dst)
{
    assert(src->capacity == dst->capacity);
    dst->len = src->len;
    copy_int_array(src->p, dst->p, src->capacity);
    copy_int_array(src->q, dst->q, src->capacity);
}

void free_sparse_array(struct sparse_array_t **S)
{
    assert(S != NULL && *S != NULL);
    free((*S)->p);
    free((*S)->q);
    free(*S);
    *S = NULL;
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
    // ??change p[j] et p[n]
    int y = S->p[n];
    S->p[n] = x;
    S->p[j] = y;
    // met q ?? jour
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


struct context_t * new_copy_context(const struct context_t *ctx)
{
    struct context_t *cpy = malloc(sizeof(*cpy));
    if (cpy == NULL)
        err(1, "impossible d'allouer un contexte");
    int nb_items = ctx->active_items->capacity;

    cpy->active_items = new_copy_sparse_array(ctx->active_items);
    cpy->active_options = malloc(nb_items * sizeof(*cpy->active_options));
    if (cpy->active_options == NULL)
        err(1, "impossible d'allouer le tableau d'options actives");
    for (int i = 0; i < nb_items; i++)
        cpy->active_options[i] = new_copy_sparse_array(ctx->active_options[i]);
    cpy->chosen_options = new_copy_int_array(ctx->chosen_options,   nb_items);
    cpy->chosen_items   = new_copy_int_array(ctx->chosen_items,     nb_items);
    cpy->child_num      = new_copy_int_array(ctx->child_num,        nb_items);
    cpy->num_children   = new_copy_int_array(ctx->num_children,     nb_items);
    cpy->level = ctx->level;
    cpy->nodes = ctx->nodes;
    cpy->solutions = ctx->solutions;

    return cpy;
}

void copy_context(const struct context_t *src, struct context_t *dst)
{
    assert(src->active_items->capacity == dst->active_items->capacity);

    int n_items = src->active_items->capacity;
    copy_sparse_array(src->active_items, dst->active_items);
    for (int i = 0; i < n_items; i++)
        copy_sparse_array(src->active_options[i], dst->active_options[i]);
    copy_int_array(src->chosen_options, dst->chosen_options,    n_items);
    copy_int_array(src->chosen_items,   dst->chosen_items,      n_items);
    copy_int_array(src->child_num,      dst->child_num,         n_items);
    copy_int_array(src->num_children,   dst->num_children,      n_items);
    dst->level = src->level;
    dst->nodes = src->nodes;
    dst->next_work_share = src->next_work_share;
    dst->solutions = src->solutions;
}

void free_context(struct context_t **ctx)
{
    assert(ctx != NULL && *ctx != NULL);
    int nb_items = (*ctx)->active_items->capacity;
    free_sparse_array(&(*ctx)->active_items);
    for (int i = 0; i < nb_items; i++)
        free_sparse_array(&(*ctx)->active_options[i]);
    free((*ctx)->chosen_options);
    free((*ctx)->child_num);
    free((*ctx)->num_children);
    free(*ctx);
    *ctx = NULL;
}

void free_instance(struct instance_t **instance)
{
   assert(instance != NULL && *instance != NULL);
   if ((*instance)->item_name != NULL)
       for (int i = 0; i < (*instance)->n_items; i++)
           free((*instance)->item_name[i]);
   free((*instance)->item_name);
   free((*instance)->options);
   free((*instance)->ptr);
   free(*instance);
   *instance = NULL;
}
