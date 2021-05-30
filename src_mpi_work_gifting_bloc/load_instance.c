#include "load_instance.h"

#include <stdio.h>
#include <err.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <mpi.h>

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
                    err(1, "erreur lors de la lecture de %s", filename);
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
                    err(1, "erreur lors de la lecture de %s", filename);
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

struct instance_t * load_instance_bcast(const char *filename)
{
    struct instance_t *instance = NULL;
    if (rank == 0) {
        instance = load_matrix(filename);
    } else {
        instance = malloc(sizeof(*instance));
        if (instance == NULL)
            err(1, "Impossible d'allouer l'instance");
    }
    MPI_Bcast(&instance->n_items, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
    MPI_Bcast(&instance->n_primary, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
    MPI_Bcast(&instance->n_options, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
    int n = instance->n_items;
    int m = instance->n_options;

    if (rank != 0) {
        instance->item_name = malloc(n * sizeof(*instance->item_name));
        instance->options = malloc(n * m * sizeof(*instance->options));
        instance->ptr = malloc((m + 1) * sizeof(*instance->ptr));

        if (instance->item_name == NULL || instance->options == NULL
                || instance->ptr == NULL)
            err(1, "Impossible d'allouer l'instance");
    }
    int *item_name_lens = malloc(n * sizeof(*item_name_lens));
    if (item_name_lens == NULL)
        err(1, "Impossible d'allouer l'instance");
    if (rank == 0) {
        for (int i = 0; i < n; i++) {
            char *item_name = instance->item_name[i];
            item_name_lens[i] = (item_name == NULL) ? 0 : strlen(item_name) + 1;
        }
    }

    MPI_Bcast(item_name_lens, n, MPI_INTEGER, 0, MPI_COMM_WORLD);
    MPI_Bcast(instance->options, n * m, MPI_INTEGER, 0, MPI_COMM_WORLD);
    MPI_Bcast(instance->ptr, m + 1, MPI_INTEGER, 0, MPI_COMM_WORLD);

    for (int i = 0; i < n; i++) {
        int len = item_name_lens[i];
        if (rank != 0) {
            char *item_name = malloc(len * sizeof(*item_name));
            if (item_name == NULL)
                err(1, "Impossible d'allouer le nom d'objet");
            instance->item_name[i] = item_name;
        }
        MPI_Bcast(instance->item_name[i], len, MPI_CHAR, 0, MPI_COMM_WORLD);
    }

    return instance;
}

struct instance_t * load_instance_ibcast(const char *filename)
{
    struct instance_t *instance = NULL;
    if (rank == 0) {
        instance = load_matrix(filename);
    } else {
        instance = malloc(sizeof(*instance));
        if (instance == NULL)
            err(1, "Impossible d'allouer l'instance");
    }
    MPI_Request n_prim_req;
    MPI_Ibcast(&instance->n_primary, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, &n_prim_req);

    MPI_Bcast(&instance->n_items, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
    MPI_Bcast(&instance->n_options, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
    int n = instance->n_items;
    int m = instance->n_options;

    if (rank != 0) {
        instance->item_name = malloc(n * sizeof(*instance->item_name));
        instance->options = malloc(n * m * sizeof(*instance->options));
        instance->ptr = malloc((m + 1) * sizeof(*instance->ptr));

        if (instance->item_name == NULL || instance->options == NULL
                || instance->ptr == NULL)
            err(1, "Impossible d'allouer l'instance");
    }
    int *item_name_lens = malloc(n * sizeof(*item_name_lens));
    if (item_name_lens == NULL)
        err(1, "Impossible d'allouer l'instance");
    if (rank == 0) {
        for (int i = 0; i < n; i++) {
            char *item_name = instance->item_name[i];
            item_name_lens[i] = (item_name == NULL) ? 0 : strlen(item_name) + 1;
        }
    }

    MPI_Request opt_req;
    MPI_Request ptr_req;
    MPI_Ibcast(instance->options, n * m, MPI_INTEGER, 0, MPI_COMM_WORLD, &opt_req);
    MPI_Ibcast(instance->ptr, m + 1, MPI_INTEGER, 0, MPI_COMM_WORLD, &ptr_req);
    MPI_Bcast(item_name_lens, n, MPI_INTEGER, 0, MPI_COMM_WORLD);

    int nb_reqs = 3 + n;
    MPI_Request *reqs = malloc(nb_reqs * sizeof(*reqs));
    if (reqs == NULL)
        err(1, "Impossible d'allouer les MPI_Request");

    for (int i = 0; i < n; i++) {
        int len = item_name_lens[i];
        if (rank != 0) {
            char *item_name = malloc(len * sizeof(*item_name));
            if (item_name == NULL)
                err(1, "Impossible d'allouer le nom d'objet");
            instance->item_name[i] = item_name;
        }
        MPI_Request *req = &reqs[3 + i];
        MPI_Ibcast(instance->item_name[i], len, MPI_CHAR, 0, MPI_COMM_WORLD, req);
    }
    reqs[0] = n_prim_req;
    reqs[1] = opt_req;
    reqs[2] = ptr_req;
    MPI_Waitall(nb_reqs, reqs, MPI_STATUSES_IGNORE);

    return instance;
}
