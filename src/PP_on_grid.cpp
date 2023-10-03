#include "prioritized_planning.hpp"

int main(int argc, char **argv) {
  int n, m, k, want = std::stoi(argv[1]);
  std::string log_filename = "";
  if (argc > 2) {
    log_filename = std::string(argv[2]);
  }
  scanf("%d%d%d", &n, &m, &k);
  std::vector<int> x0(k), y0(k), x1(k), y1(k);
  for (int i = 0; i < k; i++) {
    scanf("%d%d%d%d", &x0[i], &y0[i], &x1[i], &y1[i]);
  }
  std::vector<std::vector<int>> e(n * m);
  for (int x = 0; x < n; x++) {
    for (int y = 0; y < m; y++) {
      int id = m * x + y;
      e[id].push_back(id);
      if (0 < x) {
        e[id].push_back(id - m);
      }
      if (x + 1 < n) {
        e[id].push_back(id + m);
      }
      if (0 < y) {
        e[id].push_back(id - 1);
      }
      if (y + 1 < m) {
        e[id].push_back(id + 1);
      }
    }
  }
  PrioritizedPlanning search;
  search.edges = e;
  search.coordinates = std::vector<AnonymousBhattacharya::Point>(n * m);
  for (int x = 0; x < n; x++) {
    for (int y = 0; y < m; y++) {
      int id = m * x + y;
      search.coordinates[id].x = x;
      search.coordinates[id].y = y;
    }
  }
  search.noa = k;
  search.starts = std::vector<int>(k);
  search.goals = std::vector<int>(k);
  for (int i = 0; i < k; i++) {
    search.starts[i] = m * x0[i] + y0[i];
    search.goals[i] = m * x1[i] + y1[i];
  }
  search.want = want;
  auto plans = search.search(log_filename);
  printf("%d %d %d %d\n", n, m, k, (int)plans.size());
  for (auto &plan : plans) {
    auto &routes = plan.routes;
    auto &braid = plan.braid;
    int makespan = plan.makespan;
    printf("%d %d\n", makespan, plan.cost);
    std::vector<int> curs = search.starts;
    for (int i = 0; i < k; i++) {
      printf("%d %d ", curs[i] / m, curs[i] % m);
    }
    putchar('\n');
    for (int s = 0; s < makespan; s++) {
      for (int i = 0; i < k; i++) {
        if (routes[i].size() > s) {
          curs[i] = e[curs[i]][routes[i][s]];
        }
        printf("%d %d ", curs[i] / m, curs[i] % m);
      }
      putchar('\n');
    }
    braid.print();
  }
  return 0;
}
