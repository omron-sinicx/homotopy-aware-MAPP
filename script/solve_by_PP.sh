#!/bin/bash

for ((i=0;i<$4;i++)); do
    ../build/PP_on_grid $5 < ../build/grid_instances/$1-$2-$3-$i.txt > ../build/PP_results/$1-$2-$3-$i-$5.txt &
done
