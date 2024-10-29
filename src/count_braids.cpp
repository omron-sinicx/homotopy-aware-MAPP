#include "dynnikov.hpp"
#include <set>

int main(int argc, char **argv) {
  int width, height, number_of_agents, number_of_plans;
  scanf("%d%d%d%d", &width, &height, &number_of_agents, &number_of_plans);
  int number_of_obstacles = std::stoi(argv[1]);
  std::set<Dynnikov> braids;
  for (int T = 0; T < number_of_plans; T++) {
    int makespan;
    scanf("%d", &makespan);
    char c;
    while ((c = getchar()) != '\n')
      ;
    for (int t = 0; t <= makespan; t++) {
      for (int i = 0; i < number_of_agents; i++) {
        double x, y;
        scanf("%lf%lf", &x, &y);
      }
    }

    Dynnikov b(number_of_agents + number_of_obstacles);
    for (int t = 0; t < 2*b.N - 2; t++) {
      gmp_scanf("%Zd", &(b.cd[t]));
    }
    braids.insert(b);
  }
  printf("%d\n", (int)braids.size());
  return 0;
}
