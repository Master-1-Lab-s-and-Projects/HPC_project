#!/bin/sh

if [ $# -ne 2 ]
then
    printf "Usage: %s program result_file\n" "${0##*/}"
    exit
fi

PROGRAM="$1"
RESULTS="$2"

for f in $(cat instances.list)
do
    instance_name=$(echo $f | cut -d ',' -f 1)
    nb_solutions_expected=$(echo $f | cut -d ',' -f 2)
	echo "$instance_name" >> $RESULTS

    echo "### $instance_name ###"
	for i in $(seq 1 16)
	do
        echo "$i/16"
        res=$(mpirun -N $i "$PROGRAM" --in  "../Instances/$instance_name" --progress-report 0)
        if [ $? -ne 0 ]
        then
            echo "Error" >> $RESULTS
            continue
        fi

        nb_solutions=$(echo $res | sed 's/^.*Trouvé \([0-9]\+\) sol.*$/\1/')
        if [ $nb_solutions -ne $nb_solutions_expected -a $nb_solutions -ne 0 ]
        then
            echo "Instance: $instance_name. Expected $nb_solutions_expected solutions, found $solutions !"
            exit
        fi

        time=$(echo $res | sed 's/FINI.*en \(.*\)s/\1/')
        echo $time >> $RESULTS
	done
done
