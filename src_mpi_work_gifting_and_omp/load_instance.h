#ifndef LOAD_INSTANCE_H
#define LOAD_INSTANCE_H

#include "datastructure.h"

extern int rank;           // rang du processus (MPI)

struct instance_t * load_matrix(const char *filename);
struct instance_t * load_instance_bcast(const char *filename);
struct instance_t * load_instance_ibcast(const char *filename);

#endif
