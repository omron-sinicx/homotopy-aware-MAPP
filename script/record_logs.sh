#!/bin/bash

for ((i=0;i<$4;i++)); do
    ./PP_on_grid $5 0 logs8/$1-$2-$3-$i-$5.txt 72000000 < ../build/grid_instances/$1-$2-$3-$i.txt > /dev/null &
    ./PP_on_grid $5 1 logs7/$1-$2-$3-$i-$5.txt < ../build/grid_instances/$1-$2-$3-$i.txt > /dev/null &
done
