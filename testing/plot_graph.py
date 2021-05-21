import matplotlib.pyplot as plt
import numpy as np
import sys
import os

GRAPH_DIR="graphs/"

def plot_graph(running_times_per_prog, instance_name):

    accelerations = dict()
    efficiencies = dict()
    for prog_name,times in running_times_per_prog.items():
        fst_idx = 0 if times[0] != 0 else 1
        accelerations[prog_name] = (
            fst_idx,
            [times[fst_idx]/tps for tps in times[fst_idx:]]
        )
        efficiencies[prog_name] = (
            fst_idx,
            [acc/(i+1) for i,acc in enumerate(accelerations[prog_name][1])]
        )

    # Graphe efficacité
    fig, axs = plt.subplots(2, 1, constrained_layout=True)

    for prog_name, acc in accelerations.items():
        fst_idx, acceleration = acc
        x_coords = range(fst_idx + 1, fst_idx + len(acceleration) + 1)
        axs[0].plot(x_coords, acceleration, label=prog_name)

    axs[0].set_xlabel('Nb Proc')
    axs[0].set_ylabel('Acceleration')
    axs[0].set_title('Acceleration en fonction du nombre de processus')
    axs[0].legend()


    for prog_name, eff in efficiencies.items():
        fst_idx, efficiency = eff
        x_coords = range(fst_idx + 1, fst_idx + len(efficiency) + 1)
        axs[1].plot(x_coords, efficiency, label=prog_name)

    axs[1].set_xlabel('Nb proc')
    axs[1].set_ylabel('Efficacite')
    axs[1].set_title('Efficacite en fonction du nombre de processus')
    axs[1].legend()
    fig.suptitle(instance_name)

    plt.savefig(GRAPH_DIR + instance_name[:-2] + 'png')




if __name__ == '__main__':

    if len(sys.argv) < 3:
        print("Usage:", sys.argv[0], "nb_tests_per_instance results_file1 [results_file2 ...]")
        exit(1)

    os.mkdir(GRAPH_DIR)

    tests_per_instance = int(sys.argv[1])
    result_files = sys.argv[2:]

    results_per_instance = dict()
    for file in result_files:
        with open(file, 'r') as f:
            prog_results = dict()
            for i,l in enumerate(f.readlines()):
                if i % (tests_per_instance + 1) == 0:
                    instance = l.strip()
                    prog_results[instance] = list()
                else:
                    prog_results[instance].append(float(l.strip()))

            filename = os.path.basename(f.name)
            for instance, times in prog_results.items():
                if results_per_instance.get(instance) is None:
                    results_per_instance[instance] = dict()
                results_per_instance[instance][filename[:-len('.txt')]] = times


    for instance, times_per_prog in results_per_instance.items():
        plot_graph(times_per_prog, instance)
