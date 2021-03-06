#!/bin/sh

if [ $# -ne 2 ]
then
    printf "Usage: %s program result_file\n" "${0##*/}"
    exit
fi

PROGRAM="$1"
RESULTS="$2"

nb_proc_min=2
nb_proc_max=40

for f in $(grep -v '^#' instances.list)
do
    instance_name=$(echo $f | cut -d ',' -f 1)
    nb_solutions_expected=$(echo $f | cut -d ',' -f 2)
	echo -n "$instance_name " >> $RESULTS
	echo -n "0.0 " >> $RESULTS

    echo "### $instance_name ###"
	for i in $(seq $nb_proc_min $nb_proc_max)
	do
        echo "$i/$nb_proc_max"
        res=$(mpirun --hostfile $OAR_NODEFILE -np $i "$PROGRAM" --in  "../Instances/$instance_name" --load-method bcast --progress-report 0)
        if [ $? -ne 0 ]
        then
            echo -n "error " >> $RESULTS
            continue
        fi

        nb_solutions=$(echo $res | sed 's/^.*Trouvé \([0-9]\+\) sol.*$/\1/')
        if [ $nb_solutions -ne $nb_solutions_expected -a $nb_solutions -ne 0 ]
        then
            echo "Instance: $instance_name. Expected $nb_solutions_expected solutions, found $solutions !"
            exit
        fi

        time=$(echo $res | sed 's/FINI.*en \(.*\)s/\1/')
        echo -n "$time " >> $RESULTS
	done
    echo "" >> $RESULTS
done

if [ $(grep -vc '^#' instances.list) -ne 0 ]
then
    sed -i '$ d' $RESULTS
fi
