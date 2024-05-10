#!/bin/bash

for ((i=0;i<$4;i++)); do
    ../build/optimize_plans < ../build/PP_results3/$1-$2-$3-$i-$5.txt > ../build/PP_opt_results4/$1-$2-$3-$i-$5.txt 2> /dev/null &
    ../build/optimize_plans < ../build/PPvP_results3/$1-$2-$3-$i-$5.txt > ../build/PPvP_opt_results4/$1-$2-$3-$i-$5.txt 2> /dev/null &
    #../build/optimize_plans < ../build/CBS_results//$1-$2-$3-$i.txt > ../build/CBS_opt_results2/$1-$2-$3-$i.txt 2> /dev/null &
done
