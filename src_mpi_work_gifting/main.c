#define _POSIX_C_SOURCE 200112L
#include <unistd.h>

#include "datastructure.h"
#include "load_instance.h"
#include "utility.h"
#include "exact_cover.h"

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <stdbool.h>
#include <string.h>
#include <err.h>

#include <mpi.h>

int nb_proc = 0;
int rank = 0;

double start = 0.0;

char *in_filename = NULL;              // nom du fichier contenant la matrice
bool print_solutions = false;          // affiche chaque solution
long long report_delta = 1e6;          // affiche un rapport tous les ... noeuds
long long next_report;                 // prochain rapport affiché au noeud...
long long max_solutions = 0x7fffffffffffffff;        // stop après ... solutions

void usage(char **argv)
{
    printf("%s --in FILENAME [OPTIONS]\n\n", argv[0]);
    printf("Options:\n");
    printf("--progress-report N   display a message every N nodes (0 to disable)\n");
    printf("--print-solutions     display solutions when they are found\n");
    printf("--stop-after N        stop the search once N solutions are found\n");
    printf("--debugging           pauses the program and print the PID of the processes\n");
    printf("--load-method {bcast,ibcast,parser}  how to load the instance (default:parser)\n");
    exit(0);
}

int main(int argc, char **argv)
{
    struct option longopts[7] = {
        {"in", required_argument, NULL, 'i'},
        {"progress-report", required_argument, NULL, 'v'},
        {"print-solutions", no_argument, NULL, 'p'},
        {"stop-after", required_argument, NULL, 's'},
        {"debugging", no_argument, NULL, 'g'},
        {"load-method", required_argument, NULL, 'l'},
        {NULL, 0, NULL, 0}
    };

    bool debug_mode = false;
    const char *load_method = NULL;
    char ch;
    while ((ch = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
        switch (ch) {
            case 'i':
                in_filename = optarg;
                break;
            case 'p':
                print_solutions = true;
                break;
            case 's':
                max_solutions = atoll(optarg);
                break;
            case 'v':
                report_delta = atoll(optarg);
                break;
            case 'g':
                debug_mode = true;
                break;
            case 'l':
                load_method = optarg;
                break;
            default:
                errx(1, "Unknown option\n");
        }
    }
    if (in_filename == NULL)
        usage(argv);
    next_report = report_delta;

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nb_proc);

    struct instance_t * instance = NULL;
    if (load_method == NULL || strcmp(load_method, "parser") == 0)
        instance = load_matrix(in_filename);
    else if (strcmp(load_method, "bcast") == 0)
        instance = load_instance_bcast(in_filename);
    else if (strcmp(load_method, "ibcast") == 0)
        instance = load_instance_ibcast(in_filename);
    else
        errx(1, "loading method '%s' is invalid\n", load_method);
    struct context_t * ctx = backtracking_setup(instance);
    start = wtime();



    if (debug_mode) {
        int ifl = 0;
        char hostname[256];
        gethostname(hostname, sizeof(hostname));
        printf("PID %d on %s ready for attach, rank [%d]\n", getpid(), hostname, rank);
        fflush(stdout);
        while (ifl == 0) {
            sleep(5);
        }
    }

    start = wtime();
    if (nb_proc > 1)
        launch_parallel(instance, ctx);
    if (rank == 0)
        printf("FINI. Trouvé %lld solutions en %.3fs\n", ctx->solutions,
                wtime() - start);

    MPI_Finalize();
    exit(EXIT_SUCCESS);
}
