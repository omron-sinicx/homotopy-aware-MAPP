#!/bin/bash

mkdir -p ../build/$1_PP_opt_results
mkdir -p ../build/$1_PPvP_opt_results
mkdir -p ../build/$1_ICBS_opt_results

for ((i=0;i<$3;i++)); do
    ../build/optimize_plans ../maps/$1.map < ../build/$1_PP_results/$2-$i-$4.txt > ../build/$1_PP_opt_results/$2-$i-$4.txt 2> /dev/null &
    ../build/optimize_plans ../maps/$1.map < ../build/$1_PPvP_results/$2-$i-$4.txt > ../build/$1_PPvP_opt_results/$2-$i-$4.txt 2> /dev/null &
    ../build/optimize_plans ../maps/$1.map < ../build/$1_ICBS_results/$2-$i.txt > ../build/$1_ICBS_opt_results/$2-$i.txt 2> /dev/null &
done
