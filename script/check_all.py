for i in {0..99}
do
    echo $i
    python3 ../script/check_plan.py ../build/$1/10-$i-100.txt
done
