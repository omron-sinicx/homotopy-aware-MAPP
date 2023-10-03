#!/bin/bash

for ((i=0;i<$4;i++)); do
    python3 ../script/generate_distance_instance.py $1 $2 $3 $i > ../build/grid_instances/$1-$2-$3-$i.txt
done
