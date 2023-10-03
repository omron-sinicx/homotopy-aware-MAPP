#include "braid.hpp"
#include <set>
int main(int argc, char **argv) {
  int width, height, number_of_agents, number_of_plans;
  scanf("%d%d%d%d", &width, &height, &number_of_agents, &number_of_plans);
  std::set<Braid> braids;
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

    int l;
    scanf("%d", &l);
    Braid b;
    b.N = number_of_agents;
    for (int t = 0; t < l; t++) {
      int i;
      char s[2];
      scanf("%d %s", &i, s);
      b.add(i, s[0] == '+');
    }
    b.reduce_fully();
    braids.insert(b);
  }
  printf("%d\n", (int)braids.size());
  return 0;
}
