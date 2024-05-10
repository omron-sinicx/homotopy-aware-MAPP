#!/bin/bash

for ((i=0;i<$4;i++)); do
    ../build/PP_on_grid $5 1 < ../build/grid_instances/$1-$2-$3-$i.txt > ../build/PP_results3/$1-$2-$3-$i-$5.txt &
done
