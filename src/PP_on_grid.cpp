#include "prioritized_planning.hpp"
#include <cstring>
#include <queue>

int main(int argc, char **argv) {
  int n, m, want = std::stoi(argv[1]), braid_type = 1, time_limit = 0;
  std::string log_filename = "";
  if (argc > 2) {
    braid_type = std::stoi(argv[2]);
    if (argc > 3) {
      log_filename = std::string(argv[3]);
      if (argc > 4) {
        time_limit = std::stoi(argv[4]);
      }
    }
  }
  PrioritizedPlanning search;
  search.read_grid_map(n, m);
  int k = search.noa;
  search.want = want;
  search.braid_type = (braid_type == 0 ? generator_type : dynnikov_type);
  auto plans = search.search(log_filename, time_limit);
  if (plans[0].routes.size() < k) {
    return 0;
  }
  printf("%d %d %d %d\n", n, m, k, (int)plans.size());
  for (auto &plan : plans) {
    auto &routes = plan.routes;
    auto &braid = plan.braid;
    int makespan = plan.makespan;
    printf("%d %d\n", makespan, plan.cost);
    std::vector<int> curs = search.starts;
    for (int i = 0; i < k; i++) {
      printf("%d %d ", search.coordinates[curs[i]].x,
             search.coordinates[curs[i]].y);
    }
    putchar('\n');
    for (int s = 0; s < makespan; s++) {
      for (int i = 0; i < k; i++) {
        if (routes[i].size() > s) {
          curs[i] = search.edges[curs[i]][routes[i][s]];
        }
        printf("%d %d ", search.coordinates[curs[i]].x,
               search.coordinates[curs[i]].y);
      }
      putchar('\n');
    }
    braid->print();
  }
  return 0;
}
