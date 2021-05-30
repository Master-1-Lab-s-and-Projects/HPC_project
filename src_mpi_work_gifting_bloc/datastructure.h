#ifndef DATASTRUCTURE_H
#define DATASTRUCTURE_H

#include <stdbool.h>

struct instance_t {
    int n_items;
    int n_primary;
    int n_options;
    char **item_name;   // potentiellement NULL, sinon de taille n_items
    int *options;       // l'option i contient les objets options[ptr[i]:ptr[i+1]]
    int *ptr;           // de taille n_options + 1
};

struct sparse_array_t {
    int len;           // nombre d'éléments stockés
    int capacity;      // taille maximale
    int *p;            // contenu de l'ensemble = p[0:len]
    int *q;            // taille capacity (tout comme p)
};

struct context_t {
    struct sparse_array_t *active_items;      // objets actifs
    struct sparse_array_t **active_options;   // options actives contenant l'objet i
    int *chosen_options;                      // options choisies à ce stade
    int *chosen_items;                        // items choisis à ce stade
    int *child_num;                           // numéro du fils exploré
    int *num_children;                        // nombre de fils à explorer
    int level;                                // nombre d'options choisies
    long long nodes;                          // nombre de noeuds explorés
    long long next_work_share;                // 'date' du prochain partage de travail
    long long solutions;                      // nombre de solutions trouvées
};

static const char DIGITS[62] = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9',
    'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j',
    'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't',
    'u', 'v', 'w', 'x', 'y', 'z',
    'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J',
    'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T',
    'U', 'V', 'W', 'X', 'Y', 'Z'};

bool item_is_primary(const struct instance_t *instance, int item);
bool item_is_active(const struct context_t *ctx, int item);

struct sparse_array_t * sparse_array_init(int n);
struct sparse_array_t * new_copy_sparse_array(const struct sparse_array_t *S);
void copy_sparse_array(const struct sparse_array_t *src,
        struct sparse_array_t *dst);

bool sparse_array_membership(const struct sparse_array_t *S, int x);
bool sparse_array_empty(const struct sparse_array_t *S);
void sparse_array_add(struct sparse_array_t *S, int x);
void sparse_array_remove(struct sparse_array_t *S, int x);
void sparse_array_unremove(struct sparse_array_t *S);
void sparse_array_unadd(struct sparse_array_t *S);


struct context_t * new_copy_context(const struct context_t *ctx);
void copy_context(const struct context_t *src, struct context_t *dst);
void free_context(struct context_t **ctx);
void free_instance(struct instance_t **instance);

#endif
