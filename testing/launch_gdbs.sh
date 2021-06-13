#!/bin/sh

PID_FILE=proc_pids.list
xclip -o > $PID_FILE

for p in $(grep -o "[0-9]\{4,\}" $PID_FILE | sort)
do
    gnome-terminal --tab --title=$p -- gdb -q --command=commands.gdb --pid $p
done
