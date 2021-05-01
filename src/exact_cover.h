#ifndef EXACT_COVER_H
#define EXACT_COVER_H

#include "datastructure.h"


extern int rank;           // rang du processus (MPI)
extern int nb_proc;        // nombre de processus (MPI)

extern bool print_solutions;          // affiche chaque solution
extern long long max_solutions;        // stop après ... solutions


void launch_parallel(const struct instance_t *instance, struct context_t *ctx);

/**
 * Alloue la mémoire nécessaire aux buffers utilisés pour la communication
 * inter processus (MPI).
 */
void parallel_setup(const struct instance_t *instance);
struct context_t * backtracking_setup(const struct instance_t *instance);

#endif
