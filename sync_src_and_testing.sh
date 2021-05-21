#!/bin/bash
rsync -vr src_* nancy.g5k:hpc/
rsync -v testing/* nancy.g5k:hpc/testing/
