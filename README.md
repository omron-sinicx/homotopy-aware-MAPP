# Homotopy-Aware Multi-Agent Path Planning in Plane
This is a repository for [the following paper](https://arxiv.org/abs/2310.01945):
- Kazumi Kasaura. 2023. “Homotopy-Aware Multi-Agent Path Planning in Plane.”


## Required Library

- [yamp-cpp](https://github.com/jbeder/yaml-cpp) 0.6.0
- [Eigen3](https://eigen.tuxfamily.org/index.php) 3.3
- [g2o](https://github.com/RainerKuemmerle/g2o) 1.0.0

You can use [Dockerfile](Dockerfile) to build an environment.

## Build
```
mkdir build
cd build
cmake ..
make
```
The following commands are to be run on the build folder.

## Experiments

### Evaluation of Runtime

#### Generate instances
```
mkdir grid_instances
bash ../script/generate_grid_instances.sh 29 29 40 10
```

#### Solve
```
mkdir logs
bash ../script/record_logs.sh 29 29 40 10 100
```

#### Plot results

```
python3 ../script/show_runtimes.py
python3 ../script/show_maxcd.py
```

### Plot coordinates for random braids

```
mkdir cd_values
./dynnikov_count 3 30 100 100 > cd_values/3.txt
./dynnikov_count 5 50 100 100 > cd_values/5.txt
./dynnikov_count 10 100 100 100 > cd_values/10.txt
./dynnikov_count 100 1000 100 100 > cd_values/100.txt
./dynnikov_count 1000 10000 100 100 > cd_values/1000.txt
python3 ../script/show_cd_values.py
```

### Optimization Experiment

#### Generate instances
```
bash ../script/generate_grid_instances.sh 14 14 10 100
```

#### Solve

```
mkdir PP_results
bash ../script/solve_by_PP.sh 14 14 10 100 100
mkdir PPvP_results
bash ../script/solve_by_PPvP.sh 14 14 10 100 100
```

Although this repository does not contain any optimal MAPF solver for the Optimal One baseline,
results of such a solver must be stored in a folder named `CBS_results` for the following steps.

#### Optimization

```
mkdir PP_opt_results
mkdir PPvP_opt_results
mkdir CBS_opt_results
bash ../script/optimize_PP_and_CBS.sh 14 14 10 100 100
```

#### Plot results

```
python3 ../script/show_results.py
```

#### Show histogram
```
python3 ../script/count_braids.py > braid_counts.txt
python3 ../script/show_histogram.py
```
## License
This software is released under the MIT License, see [LICENSE](LICENSE).

## Citation
```
@misc{kasaura2023homotopyaware,
      title={Homotopy-Aware Multi-Agent Path Planning in Plane}, 
      author={Kazumi Kasaura},
      year={2023},
      eprint={2310.01945},
      archivePrefix={arXiv},
      primaryClass={cs.MA}
}
```
