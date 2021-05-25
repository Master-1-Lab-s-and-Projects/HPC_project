#!/bin/bash

sync_results() {
    rsync -v "$1.g5k:hpc/results/*" results/
}

sites=(nancy)
if [ $# -gt 0 ]; then
    sites=$@
fi

for site in $sites; do
    sync_results "$site"
done
