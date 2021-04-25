#include <ctype.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#include <err.h>
#include <getopt.h>
#include <sys/time.h>

#include <omp.h>
#include <mpi.h>

#define WAS_CHOSEN 1


int rank = 0;           // rang du processus (MPI)
int nb_proc = 0;        // nombre de processus (MPI)
int solution_tag = 1;   // solution_tag utilisé pour la communication du nb de solutions (MPI)
int proc_status_tag = 2;  // proc_status_tag utilisé pour la communication de l'avancement des proc (MPI)

double start = 0.0;

char *in_filename = NULL;              // nom du fichier contenant la matrice
bool print_solutions = false;          // affiche chaque solution
long long report_delta = 1e6;          // affiche un rapport tous les ... noeuds
long long next_report;                 // prochain rapport affiché au noeud...
long long max_solutions = 0x7fffffffffffffff;        // stop après ... solutions

int *distribution_buffer = NULL;
int *status_buffer = NULL;      // contient la len, les child_num, et les upper_bounds du context
int **procs_status = NULL;      // contient les status de tous les procs esclaves


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
    int *child_num;                           // numéro du fils exploré
    int *num_children;                        // nombre de fils à explorer
    int *lower_bounds;                        //
    int *upper_bounds;
    int level;                                // nombre d'options choisies
    long long nodes;                          // nombre de noeuds explorés
    long long solutions;                      // nombre de solutions trouvées
};

static const char DIGITS[62] = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9',
    'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j',
    'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't',
    'u', 'v', 'w', 'x', 'y', 'z',
    'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J',
    'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T',
    'U', 'V', 'W', 'X', 'Y', 'Z'};


double wtime()
{
    struct timeval ts;
    gettimeofday(&ts, NULL);
    return (double) ts.tv_sec + ts.tv_usec / 1e6;
}


void usage(char **argv)
{
    printf("%s --in FILENAME [OPTIONS]\n\n", argv[0]);
    printf("Options:\n");
    printf("--progress-report N   display a message every N nodes (0 to disable)\n");
    printf("--print-solutions     display solutions when they are found\n");
    printf("--stop-after N        stop the search once N solutions are found\n");
    exit(0);
}


bool item_is_primary(const struct instance_t *instance, int item)
{
    return item < instance->n_primary;
}


void print_option(const struct instance_t *instance, int option)
{
    if (instance->item_name == NULL)
        errx(1, "tentative d'affichage sans noms d'objet");
    for (int p = instance->ptr[option]; p < instance->ptr[option + 1]; p++) {
        int item = instance->options[p];
        printf("%s ", instance->item_name[item]);
    }
    printf("\n");
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
//#pragma omp atomic
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

void solution_found(const struct instance_t *instance, struct context_t *ctx)
{
//#pragma omp atomic
    ctx->solutions++;
    if (!print_solutions)
        return;
    printf("Trouvé une nouvelle solution au niveau %d après %lld noeuds\n",
            ctx->level, ctx->nodes);
    printf("Options : \n");
    for (int i = 0; i < ctx->level; i++) {
        int option = ctx->chosen_options[i];
        printf("+ %d : ", option);
        print_option(instance, option);
    }
    printf("\n");
    printf("----------------------------------------------------\n");
}

void cover(const struct instance_t *instance, struct context_t *ctx, int item);

void choose_option(const struct instance_t *instance, struct context_t *ctx,
        int option, int chosen_item, int was_chosen)
{
    ctx->chosen_options[ctx->level] = option;
    ctx->level++;
    for (int p = instance->ptr[option]; p < instance->ptr[option + 1]; p++) {
        int item = instance->options[p];
        if (item == chosen_item && was_chosen)
            continue;
        cover(instance, ctx, item);
    }
}

void uncover(const struct instance_t *instance, struct context_t *ctx, int item);

void unchoose_option(const struct instance_t *instance, struct context_t *ctx,
        int option, int chosen_item, int was_chosen)
{
    for (int p = instance->ptr[option + 1] - 1; p >= instance->ptr[option]; p--) {
        int item = instance->options[p];
        if (item == chosen_item && was_chosen)
            continue;
        uncover(instance, ctx, item);
    }
    ctx->level--;
}


/**
 * Returns the item associated with the active
 * option containing the least number of items
 * (the smallest active option).
 */
int choose_next_item(struct context_t *ctx)
{
    int best_item = -1;
    int best_options = 0x7fffffff;
    struct sparse_array_t *active_items = ctx->active_items;

    // TODO: small number of objects, check if to parallelize
    // active_items->len should be above ~=> 10 or so
    for (int i = 0; i < active_items->len; i++) {
        int item = active_items->p[i];
        struct sparse_array_t *active_options = ctx->active_options[item];
        int k = active_options->len;
        if (k < best_options) {
            best_item = item;
            best_options = k;
        }
    }
    return best_item;
}

void progress_report(const struct context_t *ctx)
{
    double now = wtime();
    printf("[rank %d] Exploré %lld noeuds, trouvé %lld solutions, temps écoulé %.1fs. ",
            rank, ctx->nodes, ctx->solutions, now - start);
    int i = 0;
    for (int k = 0; k < ctx->level; k++) {
        if (i > 44)
            break;
        int n = ctx->child_num[k];
        int m = ctx->num_children[k];
        if (m == 1)
            continue;
        printf("%c%c ", (n < 62) ? DIGITS[n] : '*', (m < 62) ? DIGITS[m] : '*');
        i++;
    }
    printf("\n"),
        next_report += report_delta;
}

void deactivate(const struct instance_t *instance, struct context_t *ctx,
        int option, int covered_item);

void cover(const struct instance_t *instance, struct context_t *ctx, int item)
{
    if (item_is_primary(instance, item))
        sparse_array_remove(ctx->active_items, item);
    struct sparse_array_t *active_options = ctx->active_options[item];
    for (int i = 0; i < active_options->len; i++) {
        int option = active_options->p[i];
        deactivate(instance, ctx, option, item);
    }
}


void deactivate(const struct instance_t *instance, struct context_t *ctx,
        int option, int covered_item)
{
    for (int k = instance->ptr[option]; k < instance->ptr[option+1]; k++) {
        int item = instance->options[k];
        if (item == covered_item)
            continue;
        sparse_array_remove(ctx->active_options[item], option);
    }
}


void reactivate(const struct instance_t *instance, struct context_t *ctx,
        int option, int uncovered_item);

void uncover(const struct instance_t *instance, struct context_t *ctx, int item)
{
    struct sparse_array_t *active_options = ctx->active_options[item];
    for (int i = active_options->len - 1; i >= 0; i--) {
        int option = active_options->p[i];
        reactivate(instance, ctx, option, item);
    }
    if (item_is_primary(instance, item))
        sparse_array_unremove(ctx->active_items);
}


void reactivate(const struct instance_t *instance, struct context_t *ctx,
        int option, int uncovered_item)
{
//#pragma omp parallel for
    for (int k = instance->ptr[option + 1] - 1; k >= instance->ptr[option]; k--) {
        int item = instance->options[k];
        if (item == uncovered_item)
            continue;
        sparse_array_unremove(ctx->active_options[item]);
    }
}


// TODO: extract this function in another C source file
struct instance_t * load_matrix(const char *filename)
{
    struct instance_t *instance = malloc(sizeof(*instance));
    if (instance == NULL)
        err(1, "Impossible d'allouer l'instance");
    FILE *in = fopen(filename, "r");
    if (in == NULL)
        err(1, "Impossible d'ouvrir %s en lecture", filename);
    int n_it, n_op;
    if (fscanf(in, "%d %d\n", &n_it, &n_op) != 2)
        errx(1, "Erreur de lecture de la taille du problème\n");
    if (n_it == 0 || n_op == 0)
        errx(1, "Impossible d'avoir 0 objets ou 0 options");
    instance->n_items = n_it;
    instance->n_primary = 0;
    instance->n_options = n_op;
    instance->item_name = malloc(n_it * sizeof(char *));
    instance->ptr = malloc((n_op + 1) * sizeof(int));
    instance->options = malloc(n_it * n_op *sizeof(int));         // surallocation massive
    if (instance->item_name == NULL || instance->ptr == NULL || instance->options == NULL)
        err(1, "Impossible d'allouer la mémoire pour stocker la matrice");


    enum state_t {START, ID, WHITESPACE, BAR, ENDLINE, ENDFILE};
    enum state_t state = START;

    char buffer[256];
    int i = 0;     // prochain octet disponible du buffer
    int n = 0;     // dernier octet disponible du buffer

    char id[65];
    id[64] = 0;    // sentinelle à la fin, quoi qu'il arrive
    int j = 0;     // longueur de l'identifiant en cours de lecture

    int current_item = 0;
    while (state != ENDLINE) {
        enum state_t prev_state = state;
        if (i >= n) {
            n = fread(buffer, 1, 256, in);
            if (n == 0) {
                if (feof(in)) {
                    state = ENDFILE;
                }
                if (ferror(in))
                    err(1, "erreur lors de la lecture de %s", in_filename);
            }
            i = 0;

        }
        if (state == ENDFILE) {
            // don't examine buffer[i]
        } else if (buffer[i] == '\n') {
            state = ENDLINE;
        } else if (buffer[i] == '|') {
            state = BAR;
        } else if (isspace(buffer[i])) {
            state = WHITESPACE;
        } else {
            state = ID;
        }

        // traite le caractère lu
        if (state == ID) {
            if (j == 64)
                errx(1, "nom d'objet trop long : %s", id);
            id[j] = buffer[i];
            j++;
        }
        if (prev_state == ID && state != ID) {
            id[j] = '\0';
            if (current_item == instance->n_items)
                errx(1, "Objet excedentaire : %s", id);
            for (int k = 0; k < current_item; k++)
                if (strcmp(id, instance->item_name[k]) == 0)
                    errx(1, "Nom d'objets dupliqué : %s", id);
            instance->item_name[current_item] = malloc(j+1);
            strcpy(instance->item_name[current_item], id);
            current_item++;
            j = 0;


        }
        if (state == BAR)
            instance->n_primary = current_item;
        if (state == ENDFILE)
            errx(1, "Fin de fichier prématurée");
        // passe au prochain caractère
        i++;
    }
    if (current_item != instance->n_items)
        errx(1, "Incohérence : %d objets attendus mais seulement %d fournis\n",
                instance->n_items, current_item);
    if (instance->n_primary == 0)
        instance->n_primary = instance->n_items;

    int current_option = 0;
    int p = 0;       // pointeur courant dans instance->options
    instance->ptr[0] = p;
    bool has_primary = false;
    while (state != ENDFILE) {
        enum state_t prev_state = state;
        if (i >= n) {
            n = fread(buffer, 1, 256, in);
            if (n == 0) {
                if (feof(in)) {
                    state = ENDFILE;
                }
                if (ferror(in))
                    err(1, "erreur lors de la lecture de %s", in_filename);
            }
            i = 0;

        }
        if (state == ENDFILE) {
            // don't examine buffer[i]
        } else if (buffer[i] == '\n') {
            state = ENDLINE;
        } else if (buffer[i] == '|') {
            state = BAR;
        } else if (isspace(buffer[i])) {
            state = WHITESPACE;
        } else {
            state = ID;
        }

        // traite le caractère lu
        if (state == ID) {
            if (j == 64)
                errx(1, "nom d'objet trop long : %s", id);
            id[j] = buffer[i];
            j++;
        }
        if (prev_state == ID && state != ID) {
            id[j] = '\0';
            // identifie le numéro de l'objet en question
            int item_number = -1;
            for (int k = 0; k < instance->n_items; k++)
                if (strcmp(id, instance->item_name[k]) == 0) {
                    item_number = k;
                    break;
                }
            if (item_number == -1)
                errx(1, "Objet %s inconnu dans l'option #%d", id, current_option);
            // détecte les objets répétés
            for (int k = instance->ptr[current_option]; k < p; k++)
                if (item_number == instance->options[k])
                    errx(1, "Objet %s répété dans l'option %d\n",
                            instance->item_name[item_number], current_option);
            instance->options[p] = item_number;
            p++;
            has_primary |= item_is_primary(instance, item_number);
            j = 0;


        }
        if (state == BAR) {
            errx(1, "Trouvé | dans une option.");
        }
        if ((state == ENDLINE || state == ENDFILE)) {
            // esquive les lignes vides
            if (p > instance->ptr[current_option]) {
                if (current_option == instance->n_options)
                    errx(1, "Option excédentaire");
                if (!has_primary)
                    errx(1, "Option %d sans objet primaire\n", current_option);
                current_option++;
                instance->ptr[current_option] = p;
                has_primary = false;


            }
        }
        // passe au prochain caractère
        i++;
    }
    if (current_option != instance->n_options)
        errx(1, "Incohérence : %d options attendues mais seulement %d fournies\n",
                instance->n_options, current_option);


    fclose(in);
    fprintf(stderr, "Lu %d objets (%d principaux) et %d options\n",
            instance->n_items, instance->n_primary, instance->n_options);
    return instance;
}


struct context_t * backtracking_setup(const struct instance_t *instance)
{
    struct context_t *ctx = malloc(sizeof(*ctx));
    if (ctx == NULL)
        err(1, "impossible d'allouer un contexte");
    ctx->level = 0;
    ctx->nodes = 0;
    ctx->solutions = 0;
    int n = instance->n_items;
    int m = instance->n_options;
    ctx->active_options = malloc(n * sizeof(*ctx->active_options));
    ctx->chosen_options = malloc(n * sizeof(*ctx->chosen_options));
    ctx->child_num = malloc(n * sizeof(*ctx->child_num));
    ctx->num_children = malloc(n * sizeof(*ctx->num_children));
    ctx->lower_bounds = malloc(n * sizeof(*ctx->lower_bounds));
    ctx->upper_bounds = malloc(n * sizeof(*ctx->upper_bounds));

    status_buffer = malloc((2 * n + 1) * sizeof(*status_buffer));

    if (ctx->active_options == NULL || ctx->chosen_options == NULL
            || ctx->child_num == NULL || ctx->num_children == NULL
            || ctx->lower_bounds == NULL || ctx->upper_bounds == NULL
            || status_buffer == NULL)
        err(1, "impossible d'allouer le contexte");
    ctx->active_items = sparse_array_init(n);

    // TODO: probably not possible... modification to array
    for (int item = 0; item < instance->n_primary; item++)
        sparse_array_add(ctx->active_items, item);

    // TODO: probably not possible... parallelized malloc !?
    for (int item = 0; item < n; item++)
        ctx->active_options[item] = sparse_array_init(m);

    /*
     * [o_1, ..., 0_m]
     * instance->ptr[o_i] <= k < instance->ptr[o_(i + 1)]
     * TODO: check options disjointed..
     */
    for (int option = 0; option < m; option++)
        for (int k = instance->ptr[option]; k < instance->ptr[option + 1]; k++) {
            int item = instance->options[k];
            sparse_array_add(ctx->active_options[item], option);
        }

    return ctx;
}

void print_sparse_array(struct sparse_array_t *arr) {
    printf("[");
    for (int i=0; i < arr->len; i++)
        printf("%d%s", arr->p[i], (i == arr->len - 1)? "]\n": ", ");
}

void print_context(const struct context_t *ctx) {
    printf("Context: \n");
    printf("* active_items: ");
    print_sparse_array(ctx->active_items);

    printf("* active_options: \n");
    for (int i=0; i < ctx->active_items->len; i++) {
        int item = ctx->active_items->p[i];
        printf("\t%d, ", item);
        print_sparse_array(ctx->active_options[item]);
    }
    printf("\n");
}

void build_context(const struct instance_t *instance,
        struct context_t *ctx, const int *chosen_options, int n_options)
{
    for (int i=0; i < n_options; i++) {
        choose_option(instance, ctx, chosen_options[i], -1, !WAS_CHOSEN);
    }
}

void reset_context(const struct instance_t *instance,
        struct context_t *ctx, const int *chosen_options, int n_options)
{
    for (int i=n_options - 1; i >= 0; i--) {
        unchoose_option(instance, ctx, chosen_options[i], -1, !WAS_CHOSEN);
    }
}

void collect_results(struct context_t *ctx) {
    printf("Found %lld solutions - %lld nodes - proc of rank %d\n",
            ctx->solutions, ctx->nodes, rank);
    if (rank == 0) {
        for (int i = 1; i < nb_proc; i++) {
            int rcvd_solutions;
            MPI_Recv(&rcvd_solutions, 1, MPI_INTEGER, i, solution_tag, MPI_COMM_WORLD, NULL);
            // receive nb_solution for proc of rank i
            ctx->solutions += rcvd_solutions;
        }
    } else {
        MPI_Send(&ctx->solutions, 1, MPI_INTEGER, 0, solution_tag, MPI_COMM_WORLD);
    }
}

int master_asked_status()
{
    MPI_Status status;
    MPI_Request request;
    int res;
    MPI_Irecv(&res, 1, MPI_INTEGER, 0, proc_status_tag, MPI_COMM_WORLD, &request);
    int req_completed;
    MPI_Test(&request, &req_completed, &status);

    return req_completed;
}

/**
 * Sends upper_bounds, and child_num to the master
 */
void send_status(const struct instance_t *instance, const struct context_t *ctx)
{
    int n = instance->n_items;
    status_buffer[0] = n;
    memcpy(status_buffer + 1, ctx->child_num, n * sizeof(*ctx->child_num));
    memcpy(status_buffer + 1 + n, ctx->upper_bounds, n * sizeof(*ctx->upper_bounds));
    MPI_Send(status_buffer, 2 * n + 1, MPI_INTEGER, 0, proc_status_tag, MPI_COMM_WORLD);
}

/**
 * TODO:
 * - use a different tag
 * - maybe return the number of the proc to send the data to
 */
int master_asked_distribution()
{
    MPI_Status status;
    MPI_Request request;
    int res;
    MPI_Irecv(&res, 1, MPI_INTEGER, 0, proc_status_tag, MPI_COMM_WORLD, &request);
    int req_completed;
    MPI_Test(&request, &req_completed, &status);

    return req_completed;
}

void distribute_work(const struct instance_t *instance, struct context_t *ctx)
{
    int distribution_level = 0; // get the highest level not completed
    int lower_bound = ctx->child_num[distribution_level] + 1; // First approx
    ctx->upper_bounds[distribution_level] = lower_bound;

    distribution_buffer[0] = distribution_level;
    distribution_buffer[1] = lower_bound;
    memcpy(distribution_buffer + 2, ctx->chosen_options, instance->n_items);

    MPI_Send(distribution_buffer, 2 + instance->n_items, MPI_INTEGER, 0,
            proc_status_tag, MPI_COMM_WORLD);
}

/**
 * TODO: Use omp task + threadprivate (isolate context between threads)
 *
 * lower_bound and upper_bound indicate the bounds set for the current
 * level. If -1 is specified, the usual values (0 and active_options->len)
 * are used. It is useful for the level 0 (as the lower bound is linked to
 * the rank of the process) and when using the dynamic distribution,
 * to run the solver from a certain branch of the tree.
 */
void solve(const struct instance_t *instance, struct context_t *ctx,
        int lower_bound, int upper_bound)
{
    ctx->nodes++;
    if (ctx->nodes == next_report)
        progress_report(ctx);
    if (sparse_array_empty(ctx->active_items)) {
        solution_found(instance, ctx);
        return;                         /* succès : plus d'objet actif */
    }

    if (master_asked_status()) {
        send_status(instance, ctx);
    }
    if (master_asked_distribution()) {
        distribute_work(instance, ctx);
    }


    int chosen_item = choose_next_item(ctx);
    struct sparse_array_t *active_options = ctx->active_options[chosen_item];
    if (sparse_array_empty(active_options))
        return;           /* échec : impossible de couvrir chosen_item */

    cover(instance, ctx, chosen_item);
    ctx->num_children[ctx->level] = active_options->len;


    ctx->lower_bounds[ctx->level] = (lower_bound == -1) ? 0 : lower_bound;
    ctx->upper_bounds[ctx->level] = (upper_bound == -1) ? active_options->len : upper_bound;


    for (int k = ctx->lower_bounds[ctx->level]; k < ctx->upper_bounds[ctx->level];) {
        int option = active_options->p[k];
        ctx->child_num[ctx->level] = k;
        choose_option(instance, ctx, option, chosen_item, WAS_CHOSEN);
        // On assigne les taches ici, -> idée savoir qu'il ont fini leur travail quand il sort de cette boucle avec l'id k de la ligne 564
        // /!\ On envoie au tache le contexte, après avoir choisi l'option, et on sort de la boucle quand on a fini solve au niveau k :peepoez:
        // faudra faire marcher ca

        //#pragma omp task
        solve(instance, ctx, -1, -1);

        if (ctx->solutions >= max_solutions)
            return;
        //apres on communique le nb de solutions et on attends une prochaine tacches
        // Si on veut aller plus loin -> methode arbre binomial comme ca on parallelise le plus possible mais il risque d'avoir plus de noeuds travailleurs
        unchoose_option(instance, ctx, option, chosen_item, WAS_CHOSEN);


        // On ne distribue qu'au premier niveau
        k += (ctx->level == 0) ? (nb_proc - 1) : 1;
    }

    //Typiquement on rajoute une boucle ici pour les proc 1 ... n-1 et on les fait tourner a l'infini, on les arrete avec un signal depuis le programme principal

    uncover(instance, ctx, chosen_item);                      /* backtrack */
}


int choose_proc_to_distribute(int proc_to_exclude, int nb_items) {
    int cmd = 1;
    for (int i = 1; i < nb_proc; i++) {
        if (i == proc_to_exclude)
            continue;

        MPI_Send(&cmd, 1, MPI_INTEGER, i, proc_status_tag, MPI_COMM_WORLD);
        MPI_Recv(procs_status[i - 1], 2 * nb_items + 1, MPI_INTEGER,
                i, proc_status_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    // TODO, find which process is the least 'advanced'
    return -1;
}

int wait_for_done_proc(long long *rcvd_solutions)
{
    MPI_Status status;
    MPI_Recv(rcvd_solutions, 1, MPI_LONG_LONG_INT, MPI_ANY_SOURCE,
            solution_tag, MPI_COMM_WORLD, &status);
    return status.MPI_SOURCE;
}

void distribute_work_master(const struct instance_t *instance,
        int worker, int freeloader)
{
    /**
     * To distribute work we need to get:
     * - the level 'l' to distribute
     * - the first 'l' chosen options
     * - the lower bound for the level 'l'
     */

    MPI_Recv(distribution_buffer, 2 + instance->n_items, MPI_INTEGER, worker,
            proc_status_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    MPI_Send(distribution_buffer, 2 + instance->n_items, MPI_INTEGER, freeloader,
            proc_status_tag, MPI_COMM_WORLD);
}

void launch_parallel(const struct instance_t *instance, struct context_t *ctx)
{
    if (rank == 0) {
        printf("Number of proc: %d\n", nb_proc);

        while (1) {
            long long rcvd_solutions;
            int done_proc;
            int proc_to_distribute;

            done_proc = wait_for_done_proc(&rcvd_solutions);
            ctx->solutions += rcvd_solutions;

            printf("Received %lld solutions from proc %d\n",
                    rcvd_solutions, done_proc);

            if (/* some condition on the depth */ 1) {
                proc_to_distribute = choose_proc_to_distribute(done_proc, instance->n_items);
                distribute_work_master(instance, proc_to_distribute, done_proc);
            }
        }

        printf("FINI. Trouvé %lld solutions en %.1fs\n", ctx->solutions,
                wtime() - start);

    } else {
        int lower_bound = rank - 1;
        int upper_bound = -1;
        int n_options = -1;

        distribution_buffer[0] = -1;

        // loop to re-launch when done
        do {
            if (distribution_buffer[0] != -1) {
                n_options = distribution_buffer[0];
                lower_bound = distribution_buffer[1];
                build_context(instance, ctx, distribution_buffer + 2, n_options);
            }

            solve(instance, ctx, lower_bound, upper_bound);

            MPI_Send(&ctx->solutions, 1, MPI_INTEGER, 0, solution_tag, MPI_COMM_WORLD);

            if (distribution_buffer[0] != -1) {
                reset_context(instance, ctx, distribution_buffer + 2, n_options);
            }

            MPI_Recv(distribution_buffer, 2 + instance->n_items, MPI_INTEGER, 0,
                    proc_status_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        } while(distribution_buffer[0] != -1);

    }
}

void parallel_setup(const struct instance_t *instance)
{
    int n = instance->n_items;
    if (rank == 0) {
        procs_status = malloc((nb_proc - 1) * sizeof(*procs_status));
        // TODO
        if (procs_status == NULL)
            err(1, "impossible d'allouer le buffer de status");

        for (int i=0; i < nb_proc - 1; i++) {
            procs_status[i] = malloc((2 * n + 1) * sizeof(*procs_status[i]));
            // TODO
            if (procs_status[i] == NULL)
                err(1, "impossible d'allouer le buffer de status");
        }

    } else {
        status_buffer = malloc((2 * n + 1) * sizeof(*status_buffer));
        if (status_buffer == NULL)
            err(1, "impossible d'allouer le buffer de status");
    }

    distribution_buffer = malloc((n + 2) * sizeof(*distribution_buffer));
    if (distribution_buffer == NULL)
        err(1, "impossible d'allouer le buffer des options choisies");
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

    if (rank == 0) {
        print_context(ctx);
    }

    parallel_setup(instance);
    launch_parallel(instance, ctx);

    MPI_Finalize();
    exit(EXIT_SUCCESS);
}
