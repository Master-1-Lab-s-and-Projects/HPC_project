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
	#echo "Instances/$f"
	echo "$f" >> $RESULTS
	for i in $(seq 1 16)
	do
		#echo "$i"
		mpirun -N $i "$PROGRAM" --in  "Instances/$f" --progress-report 0 | sed 's/FINI.*en \(.*\)s/\1/' >> $RESULTS
	done
done
