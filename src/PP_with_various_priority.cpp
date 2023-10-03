#include "prioritized_planning.hpp"
#include <random>

int main(int argc, char **argv) {
  int n, m, k, want = std::stoi(argv[1]),
               seed = argc > 2 ? std::stoi(argv[2]) : 0;
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
  search.want = 1;
  search.starts.resize(k);
  search.goals.resize(k);
  printf("%d %d %d %d\n", n, m, k, want);
  std::mt19937 engine(seed);
  while (want--) {
    std::vector<int> priority(k);
    std::iota(priority.begin(), priority.end(), 0);
    std::shuffle(priority.begin(), priority.end(), engine);
    for (int i = 0; i < k; i++) {
      int agent = priority[i];
      search.starts[agent] = m * x0[i] + y0[i];
      search.goals[agent] = m * x1[i] + y1[i];
    }

    auto plan = search.search()[0];

    auto &routes = plan.routes;
    auto &braid = plan.braid;
    int makespan = plan.makespan;
    printf("%d %d\n", makespan, plan.cost);
    std::vector<int> curs = search.starts;
    for (int i = 0; i < k; i++) {
      int agent = priority[i];
      printf("%d %d ", curs[agent] / m, curs[agent] % m);
    }
    putchar('\n');
    for (int s = 0; s < makespan; s++) {
      for (int i = 0; i < k; i++) {
        int agent = priority[i];
        if (routes[agent].size() > s) {
          curs[agent] = e[curs[agent]][routes[agent][s]];
        }
        printf("%d %d ", curs[agent] / m, curs[agent] % m);
      }
      putchar('\n');
    }
    braid.print();
  }
  return 0;
}
