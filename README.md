# Homotopy-Aware Multi-Agent Path Planning in Plane
This is a repository for [the following paper](https://arxiv.org/abs/2310.01945):
- Kazumi Kasaura. 2023. “Homotopy-Aware Multi-Agent Path Planning in Plane.”


## Required Library

- [yamp-cpp](https://github.com/jbeder/yaml-cpp) 0.8.0
- [Eigen3](https://eigen.tuxfamily.org/index.php) 3.3
- [g2o](https://github.com/RainerKuemmerle/g2o) 1.0.0
- [gmp](https://gmplib.org/) 6.1.2

You can use [Dockerfile](Dockerfile) to build an environment.

## Build
```
mkdir build
cd build
cmake ..
make
```
The following commands are to be run on the build folder.

```
./PP_on_grid H < input_file
```
where `H` is the number of required solutions.

`sample` folder contains samples of input files.

## Experiments

The scripts for experiments are stored in `scripts` folder
The scripts to show graphs are also contained in `scripts` folder.
All scripts must be run in the build folder.

### Evaluation of Runtime

#### Generate instances
```
bash ../script/generate_map_instances.sh empty-48-48 500 10
bash ../script/generate_map_instances.sh den312d 500 10
bash ../script/generate_map_instances.sh random-64-64-10 500 10
```

#### Solve
```
bash ../script/record_logs_map.sh empty-48-48 500 10 100
bash ../script/record_logs_map.sh den312d 500 10 100
bash ../script/record_logs_map.sh random-64-64-10 500 10 100
```

#### Plot results

```
python3 ../script/show_runtimes_map.py
python3 ../script/show_maxcd.py
```

### Plot coordinates for random braids

```
bash ../script/run_cd_values.sh
python3 ../script/show_cd_values.py
```

### Optimization Experiment

#### Generate instances
```
bash ../script/generate_map_instances_2.sh empty-14-14 10 100
bash ../script/generate_map_instances_2.sh obstacle-14-14 10 100
```

#### Solve

```
bash ../script/solve_by_PP.sh empty-14-14 10 100 100
bash ../script/solve_by_PPvP.sh empty-14-14 10 100 100
bash ../script/solve_by_PP.sh obstacle-14-14 10 100 100
bash ../script/solve_by_PPvP.sh obstacle-14-14 10 100 100
```

Although this repository does not contain any optimal MAPF solver for the Optimal One baseline,
results of such a solver must be stored in folders named `empty-14-14_ICBS_results` and `obstacle-14-14_ICBS_results` for the following steps.

#### Optimization

```
bash ../script/optimize_PP_and_CBS.sh empty-14-14 10 100 100
bash ../script/optimize_PP_and_CBS.sh obstacle-14-14 10 100 100
```

#### Plot results

```
python3 ../script/show_results.py
python3 ../script/show_winning_rate.py
```

#### Show histogram
```
python3 ../script/count_braids.py empty-14-14 0 > count_empty-14-14.txt
python3 ../script/count_braids.py obstacle-14-14 5 > count_obstacle-14-14.txt
python3 ../script/show_histogram.py
```

## License
This software is released under the MIT License, see [LICENSE](LICENSE).

## Citation
```
@misc{kasaura2024homotopyaware,
      title={Homotopy-Aware Multi-Agent Path Planning in Plane}, 
      author={Kazumi Kasaura},
      year={2024},
      eprint={2310.01945},
      archivePrefix={arXiv},
      primaryClass={cs.MA}
}
```
