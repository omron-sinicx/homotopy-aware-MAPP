# Required Library

- [yamp-cpp](https://github.com/jbeder/yaml-cpp) 0.6.0
- [Eigen3](https://eigen.tuxfamily.org/index.php) 3.3
- [g2o](https://github.com/RainerKuemmerle/g2o) 1.0.0

You can use [Dockerfile](Dockerfile) to build an environment.

# Build
```
mkdir build
cd build
cmake ..
make
```
The following commands are to be run on the build folder.

# Experiments

## Evaluation of Runtime

### Generate instances
```
mkdir grid_instances
bash ../script/generate_grid_instances.sh 29 29 40 10
```

### Solve
```
mkdir logs
bash ../script/record_logs.sh 29 29 40 10 100
```

### Plot results

```
python3 ../script/show_makespan.py
python3 ../script/show_lengths.py
python3 ../script/show_runtimes.py
```

## Optimization Experiment

### Generate instances
```
bash ../script/generate_grid_instances.sh 14 14 10 100
```

### Solve

```
mkdir PP_results
bash ../script/solve_by_PP.sh 14 14 10 100 100
mkdir PPvP_results
bash ../script/solve_by_PPvP.sh 14 14 10 100 100
```

Although this repository does not contain any optimal MAPF solver for the Optimal One baseline,
results of such a solver must be stored in a folder named `CBS_results` for the following steps.

### Optimization

```
mkdir PP_opt_results
mkdir PPvP_opt_results
mkdir CBS_opt_results
bash ../script/optimize_PP_and_CBS.sh 14 14 10 100 100
```

### Plot results

```
python3 ../script/show_results.py
```

### Show histogram
```
python3 ../script/count_braids.py > braid_counts.txt
python3 ../script/show_histogram.py
```
