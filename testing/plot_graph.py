#!/usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np
import sys
import os

GRAPH_DIR="graphs/"

def plot_graph(running_times_per_prog, instance_name):

    accelerations = dict()
    efficiencies = dict()
    times_per_prog = dict()
    for prog_name, times in running_times_per_prog.items():
        idxs = [i for i,t in enumerate(times) if t != 0 and t != "error"]
        times = [times[i] for i in idxs]
        accelerations[prog_name] = (
            idxs,
            [times[0]/tps for tps in times]
        )
        efficiencies[prog_name] = (
            idxs,
            [acc/(i+1 - idxs[0]) for i,acc in zip(*accelerations[prog_name])]
        )
        times_per_prog[prog_name] = (
            idxs,
            times
        )

    fig, axs = plt.subplots(2, 1, constrained_layout=True)

    # Subgraph acceleration
    for prog_name, acc in accelerations.items():
        axs[0].plot(acc[0], acc[1], label=prog_name)

    all_idxs = [idxs for idxs,_ in accelerations.values()]
    max_nb_procs = max(map(max, all_idxs)) + 1
    x_coords = [*range(1, max_nb_procs + 1)]
    axs[0].plot(x_coords)

    axs[0].set_xlabel('Nb threads')
    axs[0].set_ylabel('Acceleration')
    axs[0].set_title('Acceleration en fonction du nombre de threads')
    axs[0].grid(True, axis='y', linestyle='--')
    axs[0].legend()


    # Subgraph efficiency
    for prog_name, eff in efficiencies.items():
        axs[1].plot(eff[0], eff[1], label=prog_name)

    axs[1].set_xlabel('Nb threads')
    axs[1].set_ylabel('Efficacite')
    axs[1].set_title('Efficacite en fonction du nombre de threads')
    axs[1].grid(True, axis='y', linestyle='--')
    axs[1].legend()

    fig.suptitle(instance_name)
    plt.savefig(GRAPH_DIR + instance_name[:-2] + 'png')


    # Graph Times
    fig, axs = plt.subplots(1, 1, constrained_layout=True)
    for prog_name,times in times_per_prog.items():
        axs.plot(times[0], times[1], label=prog_name)

    axs.set_xlabel('Nb Proc')
    axs.set_ylabel('Temps (s)')
    axs.set_title("Temps d'exécution en fonction du nombre de processus")
    axs.grid(True, axis='y', linestyle='--')
    axs.legend()

    fig.suptitle(instance_name)
    plt.savefig(GRAPH_DIR + 'temps_' + instance_name[:-2] + 'png')





if __name__ == '__main__':

    if len(sys.argv) < 3:
        print("Usage:", sys.argv[0], "[-d output_directory] results_file1 [results_file2 ...]")
        exit(1)

    args = sys.argv[1:]
    option_pos = next((i for i,x in enumerate(args) if x == '-d'), None)
    if option_pos is not None:
        GRAPH_DIR = args[option_pos + 1]
        if GRAPH_DIR[-1] != '/':
            GRAPH_DIR += '/'
        args = args[:option_pos] + args[option_pos + 2:]

    if not os.path.isdir(GRAPH_DIR):
        os.mkdir(GRAPH_DIR)

    result_files = args

    results_per_instance = dict()
    for file in result_files:
        with open(file, 'r') as f:
            prog_results = dict()
            for l in f.readlines():
                line = l.split()
                instance = line[0]
                prog_results[instance] = [*map(lambda x: x if x == "error" else float(x), line[1:])]

            filename = os.path.basename(f.name)
            for instance, times in prog_results.items():
                if results_per_instance.get(instance) is None:
                    results_per_instance[instance] = dict()
                results_per_instance[instance][filename[:-len('.txt')]] = times


    for instance, times_per_prog in results_per_instance.items():
        plot_graph(times_per_prog, instance)
