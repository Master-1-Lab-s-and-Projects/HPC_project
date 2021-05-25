#define _POSIX_C_SOURCE 200112L
#include <unistd.h>

#include "datastructure.h"
#include "matrix_parser.h"
#include "utility.h"
#include "exact_cover.h"

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <stdbool.h>
#include <err.h>

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
    bool debug_mode = false;
    struct option longopts[6] = {
        {"in", required_argument, NULL, 'i'},
        {"progress-report", required_argument, NULL, 'v'},
        {"print-solutions", no_argument, NULL, 'p'},
        {"stop-after", required_argument, NULL, 's'},
        {"debugging", no_argument, NULL, 'g'},
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
            case 'g':
                debug_mode = true;
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












    if (debug_mode) {
        int ifl = 0;
        char hostname[256];
        gethostname(hostname, sizeof(hostname));
        printf("PID %d on %s ready for attach\n", getpid(), hostname);
        fflush(stdout);
        while (ifl == 0) {
            sleep(5);
        }
    }

    start = wtime();
    launch_parallel(instance, ctx);
    printf("FINI. Trouvé %lld solutions en %.1fs\n", ctx->solutions,
            wtime() - start);

    free_instance(&instance);
    free_context(&ctx);

    exit(EXIT_SUCCESS);
}
