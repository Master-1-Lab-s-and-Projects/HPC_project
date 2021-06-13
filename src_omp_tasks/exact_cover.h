#ifndef EXACT_COVER_H
#define EXACT_COVER_H

#include "datastructure.h"

extern bool print_solutions;          // affiche chaque solution
extern long long max_solutions;        // stop après ... solutions


void launch_parallel(const struct instance_t *instance, struct context_t *ctx);
struct context_t * backtracking_setup(const struct instance_t *instance);

#endif
