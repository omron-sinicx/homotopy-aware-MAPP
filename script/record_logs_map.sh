#!/bin/bash

for ((i=0;i<$3;i++)); do
    ./PP_on_grid $4 0 logs_$1_deh/$2-$i-$4.txt 72000000 < ../build/$1_instances/$2-$i.txt > /dev/null &
    ./PP_on_grid $4 1 logs_$1/$2-$i-$4.txt < ../build/$1_instances/$2-$i.txt > /dev/null &
done
