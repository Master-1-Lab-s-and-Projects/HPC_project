#!/bin/sh

INSTANCE_LIST=instances.list

for instance in $(cat $INSTANCE_LIST)
do
    echo "###### $instance ######"
    mpirun -N 6 --oversubscribe ../src/exact_cover.exe --in ../Instances/$instance  --progress-report 0 | grep -v 'UCX  WARN  object'
    if [ $? -ne 0 ]
    then
        break
    fi

done
