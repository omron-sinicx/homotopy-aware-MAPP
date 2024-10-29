#!/bin/bash

mkdir -p ../build/$1_PPvP_results
for ((i=0;i<$3;i++)); do
    ../build/PP_with_various_priority $4 $i < ../build/$1_instances/$2-$i.txt > ../build/$1_PPvP_results/$2-$i-$4.txt &
done
