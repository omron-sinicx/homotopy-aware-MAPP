
#!/bin/bash

for ((i=0;i<$4;i++)); do
    ../build/PP_on_grid $5 logs/$1-$2-$3-$i-$5.txt < ../build/grid_instances/$1-$2-$3-$i.txt > /dev/null
done
