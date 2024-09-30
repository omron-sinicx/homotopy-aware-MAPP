#!/bin/bash

for ((i=0;i<$3;i++)); do
    python3 ../script/generate_grid_map_instance.py ../maps/$1.map $2 $i > ../build/$1_instances/$2-$i.txt &
done
