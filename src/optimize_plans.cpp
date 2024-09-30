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

void input_grid(
    const std::string &file_path,
    std::shared_ptr<std::vector<std::vector<Eigen::Vector2d>>> &boundaries) {
  FILE *in = fopen(file_path.c_str(), "r");
  if (in == NULL) {
    throw std::runtime_error("Map file cannot be opened\n");
  }
  int height, width;
  fscanf(in, "type octile\nheight %d\nwidth %d\nmap\n", &height, &width);
  std::vector<std::vector<bool>> grid(
      height + 2,
      std::vector<bool>(width + 2,
                        false)); // fill boundaries by unpassable blocks
  for (int i = 0; i < height; i++) {
    char c;
    for (int j = 0; j < width; j++) {
      c = getc(in);
      grid[height - i][j + 1] = c == '.' || c == 'G';
    }
    while ((c = getc(in)) == '\n' || c == '\r')
      ;
    ungetc(c, in);
  }
  fclose(in);

  std::vector<std::vector<bool>> visited(height + 2,
                                         std::vector<bool>(width + 2, false));
  bool found_outer = false;
  const int dx[4] = {1, 0, -1, 0},
            dy[4] = {0, 1, 0, -1}; // left:0 up:1 right:2 dowm:3
  const int sx[4] = {0, -1, -1, 0}, sy[4] = {0, 0, -1, -1};
  for (int start_y = 1; start_y <= height; start_y++) {
    for (int start_x = 1; start_x <= width; start_x++) {
      if (!(!grid[start_y][start_x - 1] && grid[start_y][start_x]) ||
          visited[start_y][start_x]) {
        continue;
      }
      // left hand method
      std::vector<int> xs, ys;
      int x = start_x, y = start_y;
      int direction = 1;
      do {
        if (direction == 1)
          visited[y][x] = true;
        x += dx[direction];
        y += dy[direction];
        int ldir = (direction + 3) % 4, rdir = (direction + 1) % 4;
        if (!grid[y + sy[ldir]][x + sx[ldir]]) {
          xs.push_back(x);
          ys.push_back(y);
          direction = ldir;
        } else if (grid[y + sy[direction]][x + sx[direction]]) {
          xs.push_back(x);
          ys.push_back(y);
          direction = rdir;
        }
      } while (x != start_x || y != start_y);

      std::vector<Eigen::Vector2d> obstacle(xs.size());
      for (int i = 0; i < xs.size(); i++) {
        obstacle[i] = Eigen::Vector2d(xs[i] - 1.5, ys[i] - 1.5);
      }
      boundaries->push_back(obstacle);
    }
  }
}

int main(int argc, char **argv) {
  int width, height, number_of_agents, number_of_plans;
  scanf("%d%d%d%d", &width, &height, &number_of_agents, &number_of_plans);
  std::string map_path(argv[1]);
  std::string config_file_path(argc > 2 ? argv[2]
                                        : "../config/optimizer_config.yaml");
  auto config = YAML::LoadFile(config_file_path);
  std::shared_ptr<std::vector<std::vector<Eigen::Vector2d>>> boundaries(
      new std::vector<std::vector<Eigen::Vector2d>>);
  input_grid(map_path, boundaries);
  TrajectoryOptimizer optimizer(config, number_of_agents, boundaries);
  bool output_plan = config["output_plan"].as<bool>();
  if (output_plan) {
    printf("%d %d %d %d\n", width, height, number_of_agents, number_of_plans);
  } else {
    printf("%d %d\n", number_of_agents, number_of_plans);
  }
  // std::default_random_engine engine(0);
  // double move_distance =
  //    config["move_distance"] ? config["move_distance"].as<double>() : 0.5;
  // std::uniform_real_distribution<> dist(-move_distance, move_distance);
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
