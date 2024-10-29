#!/bin/bash

mkdir -p logs_HR_$1
mkdir -p logs_Dyn_$1
for ((i=0;i<$3;i++)); do
    ./PP_on_grid $4 0 logs_HR_$1/$2-$i-$4.txt 72000000 < ../build/$1_instances/$2-$i.txt > /dev/null &
    ./PP_on_grid $4 1 logs_Dyn_$1/$2-$i-$4.txt < ../build/$1_instances/$2-$i.txt > /dev/null &
done
