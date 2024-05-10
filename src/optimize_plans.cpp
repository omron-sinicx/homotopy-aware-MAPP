#include "optimizer.cpp"
#include <random>

#include "braid.hpp"

Braid calculate_braid(const std::vector<std::vector<Eigen::Vector2d>> &plan) {
  int number_of_agents = plan.size(), number_of_way_points = plan[0].size() - 1;
  Braid b;
  b.N = number_of_agents;
  std::vector<std::pair<double, int>> start_indexes(number_of_agents);
  for (int i = 0; i < number_of_agents; i++) {
    start_indexes[i] = std::make_pair(plan[i][0].x() + EPS * plan[i][0].y(), i);
  }
  std::sort(start_indexes.begin(), start_indexes.end());
  std::vector<int> ids(number_of_agents);
  for (int i = 0; i < number_of_agents; i++) {
    ids[i] = start_indexes[i].second;
  }
  for (int t = 1; t <= number_of_way_points; t++) {
    for (int i = 0; i + 1 < number_of_agents; i++) {
      if (plan[ids[i]][t].x() + EPS * plan[ids[i]][t].y() >
          plan[ids[i + 1]][t].x() + EPS * plan[ids[i + 1]][t].y()) {
        b.add(i, plan[ids[i]][t].y() > plan[ids[i + 1]][t].y());
        std::swap(ids[i], ids[i + 1]);
      }
    }
  }
  return b;
}

int main(int argc, char **argv) {
  int width, height, number_of_agents, number_of_plans;
  scanf("%d%d%d%d", &width, &height, &number_of_agents, &number_of_plans);
  std::string config_file_path(argc > 1 ? argv[1]
                                        : "../config/optimizer_config.yaml");
  auto config = YAML::LoadFile(config_file_path);
  TrajectoryOptimizer optimizer(config, number_of_agents);
  bool output_plan = config["output_plan"].as<bool>();
  if (output_plan) {
    printf("%d %d %d %d\n", width, height, number_of_agents,
           number_of_plans);
  } else {
    printf("%d %d\n", number_of_agents,
           number_of_plans);
  }
  std::default_random_engine engine(0);
  double move_distance =
      config["move_distance"] ? config["move_distance"].as<double>() : 0.5;
  std::uniform_real_distribution<> dist(-move_distance, move_distance);
  for (int T = 0; T < number_of_plans; T++) {
    int makespan;
    scanf("%d", &makespan);
    char c;
    while ((c = getchar()) != '\n')
      ;
    std::vector<std::vector<Eigen::Vector2d>> plan(
        number_of_agents, std::vector<Eigen::Vector2d>(makespan + 1));
    for (int t = 0; t <= makespan; t++) {
      for (int i = 0; i < number_of_agents; i++) {
        double x, y;
        scanf("%lf%lf", &x, &y);
        plan[i][t] = Eigen::Vector2d((double)x, (double)y);
      }
    }

    std::vector<Eigen::Vector2d> current_starts(number_of_agents),
        current_goals(number_of_agents);
    for (int i = 0; i < number_of_agents; i++) {
      current_starts[i] = plan[i][0];
      current_goals[i] = plan[i][makespan];
    }

    if (!output_plan && T == 0) {
      for (int i = 0; i < number_of_agents; i++) {
        printf("%lf %lf ", plan[i][0].x(), plan[i][0].y());
      }
      putchar('\n');
      for (int i = 0; i < number_of_agents; i++) {
        printf("%lf %lf ", plan[i][makespan].x(), plan[i][makespan].y());
      }
      putchar('\n');
    }
    optimizer.optimize(plan);
    auto optimized_plan = optimizer.getPlan();
    double delta = optimizer.getDelta(), cost = optimizer.get_cost();
    fprintf(stderr, "%d %lf\n", T, cost);
    if (output_plan) {
      printf("%d %lf %lf\n", optimizer.number_of_way_points, delta, cost);
      for (int t = 0; t <= optimizer.number_of_way_points; t++) {
        for (int i = 0; i < number_of_agents; i++) {
          printf("%lf %lf ", optimized_plan[i][t].x(),
                 optimized_plan[i][t].y());
        }
        putchar('\n');
      }
    }
    while ((c = getchar()) != '\n')
      ;
    while ((c = getchar()) != '\n')
      putchar(c);
    putchar('\n');
    if (!output_plan) {
      printf("%lf\n", cost);
    }
  }
  return 0;
}
