#!/bin/bash

sync_src_and_testing() {
    rsync -vr src_omp_tasks  "$1.g5k:hpc/"
}

sites=(nancy)
if [ $# -gt 0 ]; then
    sites=$@
fi

for site in $sites; do
    sync_src_and_testing "$site"
done
