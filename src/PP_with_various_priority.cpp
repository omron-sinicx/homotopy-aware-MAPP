#include "prioritized_planning.hpp"
#include <random>

int main(int argc, char **argv) {
  int n, m, k, want = std::stoi(argv[1]),
               seed = argc > 2 ? std::stoi(argv[2]) : 0;
  PrioritizedPlanning search;
  search.read_grid_map(n, m);
  k = search.noa;
  search.want = 1;
  auto original_starts = search.starts, original_goals = search.goals;
  printf("%d %d %d %d\n", n, m, k, want);
  std::mt19937 engine(seed);
  while (want) {
    std::vector<int> priority(k);
    std::iota(priority.begin(), priority.end(), 0);
    std::shuffle(priority.begin(), priority.end(), engine);
    for (int i = 0; i < k; i++) {
      search.starts[priority[i]] = original_starts[i];
      search.goals[priority[i]] = original_goals[i];
    }

    PrioritizedPlanning::Plan plan;

    try {
      plan = search.search()[0];
    } catch (const std::runtime_error &e) {
      continue;
    }
    want--;

    auto &routes = plan.routes;
    auto &braid = plan.braid;
    int makespan = plan.makespan;
    printf("%d %d\n", makespan, plan.cost);
    std::vector<int> curs = search.starts;
    for (int i = 0; i < k; i++) {
      int agent = priority[i];
      printf("%d %d ", search.coordinates[curs[agent]].x,
             search.coordinates[curs[agent]].y);
    }
    putchar('\n');
    for (int s = 0; s < makespan; s++) {
      for (int i = 0; i < k; i++) {
        int agent = priority[i];
        if (routes[agent].size() > s) {
          curs[agent] = search.edges[curs[agent]][routes[agent][s]];
        }
        printf("%d %d ", search.coordinates[curs[agent]].x,
               search.coordinates[curs[agent]].y);
      }
      putchar('\n');
    }
    braid->print();
  }
  return 0;
}
