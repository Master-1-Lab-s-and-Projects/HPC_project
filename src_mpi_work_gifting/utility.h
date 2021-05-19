#ifndef UTILITY_H
#define UTILITY_H

#include "datastructure.h"

#ifdef DEBUG
    #define DPRINTF(...) printf(__VA_ARGS__)
#else
    #define DPRINTF(...)
#endif

#define min(a,b) \
    ({ __typeof__ (a) _a = (a); \
     __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

extern int rank;                // rang du processus (MPI)
extern double start;            // date de lancement du processus
extern long long report_delta;  // affiche un rapport tous les ... noeuds
extern long long next_report;   // prochain rapport affiché au noeud...

double wtime();

void print_option(const struct instance_t *instance, int option);
void print_sparse_array(struct sparse_array_t *arr);
void print_context(const struct context_t *ctx);

void progress_report(const struct context_t *ctx);
#endif
