#!/bin/bash

mkdir -p ../build/$1_PP_results
for ((i=0;i<$3;i++)); do
    ../build/PP_on_grid $4 1 < ../build/$1_instances/$2-$i.txt > ../build/$1_PP_results/$2-$i-$4.txt &
done
