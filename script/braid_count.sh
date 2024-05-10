#!/bin/bash

Ns=(3 5 10 100 1000)

for i in ${Ns[@]}; do
    ../build/dynnikov_count $i $((10*i)) 100 100 > ../build/cd_values/$i.txt &
done
