#include "datastructure.h"
#include "matrix_parser.h"
#include "utility.h"
#include "exact_cover.h"

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <stdbool.h>
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
    exit(0);
}

int main(int argc, char **argv)
{
    struct option longopts[5] = {
        {"in", required_argument, NULL, 'i'},
        {"progress-report", required_argument, NULL, 'v'},
        {"print-solutions", no_argument, NULL, 'p'},
        {"stop-after", required_argument, NULL, 's'},
        {NULL, 0, NULL, 0}
    };
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
            default:
                errx(1, "Unknown option\n");
        }
    }
    if (in_filename == NULL)
        usage(argv);
    next_report = report_delta;


    struct instance_t * instance = load_matrix(in_filename);
    struct context_t * ctx = backtracking_setup(instance);
    start = wtime();

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nb_proc);

    parallel_setup(instance);
    launch_parallel(instance, ctx);

    MPI_Finalize();
    exit(EXIT_SUCCESS);
}
