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
            [(times[fst_idx]/tps)*100 for tps in times[fst_idx:]]
        )
        efficiencies[prog_name] = (
            fst_idx,
            [(acc/(i+1))/100 for i,acc in enumerate(accelerations[prog_name][1])]
        )

    # Graphe efficacité
    fig, axs = plt.subplots(2, 1, constrained_layout=True)

    for prog_name, acc in accelerations.items():
        fst_idx, acceleration = acc
        x_coords = range(fst_idx + 1, fst_idx + len(acceleration) + 1)
        axs[0].plot(x_coords, acceleration, label=prog_name)

    axs[0].set_xlabel('Nb Proc')
    axs[0].set_ylabel('% Acceleration')
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
    d = {'bell12.ec': ['2.6', '2.4', '2.7', '2.4', '2.5', '2.4', '2.8', '2.7', '2.5', '2.5', '2.5', '2.9', '2.7', '2.6', '2.7', '2.7'],
         'bell13.ec': ['15.3', '14.8', '14.3', '14.3', '15.0', '14.2', '15.5', '15.9', '14.6', '16.0', '14.9', '15.6', '15.2', '15.5', '14.8', '15.1'],
         'bell14.ec': ['100.9', '103.0', '102.1', '104.6', '122.4', '106.9', '104.1', '113.7', '121.7', '105.7', '108.6', '104.3', '106.0', '107.0', '127.7', '107.0'],
         'matching10.ec': ['172.7', '92.1', '60.6', '65.8', '63.7', '37.5', '66.3', '50.4', '29.8', '97.8', '87.4', '78.3', '70.3', '65.3', '51.5', '40.3'],
         'matching8.ec': ['1.0', '0.7', '0.7', '0.8', '0.6', '0.7', '0.6', '0.9', '0.8', '0.9', '0.8', '0.7', '0.7', '0.6', '0.7', '1.2'],
         'matching9.ec': ['9.7', '5.0', '4.1', '3.2', '3.1', '4.1', '3.4', '2.3', '5.7', '5.1', '4.6', '3.9', '3.6', '2.9', '2.4', '1.8'],
         'pentomino_6_10.ec': ['16.2', '8.3', '6.3', '4.9', '4.4', '3.6', '4.6', '3.2', '4.5', '2.9', '6.9', '5.5', '5.0', '3.7', '2.6', '2.0'],
         'pento_plus_tetra_2x4x10.ec': ['1003.0', '519.8', '370.0', '266.6', '263.9', '253.5', '279.4', '150.1', '283.8', '170.1', '425.7', '340.3', '302.4', '216.3', '149.7', '80.3']}

    d2 = {'bell12.ec': ['2.6', '1.7', '1.4', '1.2', '1.2', '1.2', '1.1', '1.1', '1.2', '1.1', '1.1', '1.1', '1.2', '1.2', '1.2', '1.4'],
          'bell13.ec': ['15.3', '9.1', '6.9', '6.2', '5.4', '5.7', '4.7', '4.7', '4.6', '4.3', '4.7', '4.3', '4.5', '4.0', '4.0', '4.0'],
          'bell14.ec': ['122.5', '60.5', '45.9', '40.0', '35.5', '32.8', '30.8', '34.2', '27.8', '30.7', '27.6', '26.4', '26.0', '24.2', '25.1', '24.9'],
          'matching10.ec': ['162.9', '87.8', '60.2', '45.5', '36.7', '36.3', '29.8', '32.6', '29.4', '21.7', '21.7', '20.1', '21.4', '22.6', '23.1', '21.7'],
          'matching8.ec': ['1.0', '0.7', '0.7', '0.7', '0.6', '0.6', '0.7', '0.6', '0.6', '0.7', '0.7', '0.7', '0.6', '0.7', '0.6', '0.7'],
          'matching9.ec': ['8.8', '5.1', '3.5', '2.9', '2.7', '2.1', '2.3', '2.2', '2.0', '1.7', '1.7', '1.7', '1.8', '1.8', '1.8', '1.8'],
          'pentomino_6_10.ec': ['16.9', '8.3', '6.2', '4.8', '4.3', '3.6', '3.9', '3.0', '3.3', '2.8', '2.7', '2.6', '3.1', '3.1', '2.2', '2.5'],
          'pento_plus_tetra_2x4x10.ec': ['1023.0', '524.7', '355.2', '423.0', '251.0', '207.3', '194.6', '285.4', '191.3', '154.0', '171.0', '196.5', '181.5', '125.2', '125.1', '154.2']}

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
