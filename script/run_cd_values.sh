#!/bin/bash

mkdir -p ../build/cd_values

for n in 3 5 10 100 1000
do
    ./dynnikov_count $n $(($n*10)) 100 100 > ../build/cd_values/$n.txt
done
