#!/bin/bash

for ((i=0;i<$4;i++)); do
    ../build/PP_with_various_priority $5 $i < ../build/grid_instances/$1-$2-$3-$i.txt > ../build/PPvP_results3/$1-$2-$3-$i-$5.txt &
done
