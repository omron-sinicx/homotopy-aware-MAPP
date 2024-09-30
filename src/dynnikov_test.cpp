#include "braid.hpp"
#include "dynnikov.hpp"
#include <algorithm>
#include <cassert>

void test() {
  int n, N, l, c;
  scanf("%d%d%d%d", &n, &N, &l, &c);
  std::vector<Braid> braids(N, Braid(n));
  std::vector<Dynnikov> dyns(N, Dynnikov(n));
  for (int t = 0; t < N; t++) {
    for (int j = 0; j < l; j++) {
      int i = rand() % (n - 1);
      bool sign = rand() % 2;
      if (c % 2) {
        braids[t].add(i, sign);
      }
      if (c / 2 % 2) {
        dyns[t].add(i, sign);
      }
    }
    if (c % 2) {
      braids[t].reduce_fully();
    }
  }
  int count = 0;
  for (int t = 0; t < N; t++) {
    for (int s = t + 1; s < N; s++) {
      bool res;
      if (c % 2) {
        res = (braids[t] < braids[s]);
      }
      if (c / 2 % 2) {
        res = (dyns[t] < dyns[s]);
      }
      if (c % 2 && c / 2 % 2) {
        bool eq1 = !(braids[t] < braids[s]) && !(braids[s] < braids[t]);
        bool eq2 = dyns[t] == dyns[s];
        if (eq1 != eq2) {
          printf("%d %d\n", eq1, eq2);
          braids[t].print();
          braids[s].print();
          dyns[t].print();
          dyns[s].print();
        }
        assert(eq1 == eq2);
        count += eq1;
      }
    }
  }
  printf("%d\n", count);
}

void convert() {
  int n, m;
  scanf("%d%d", &n, &m);
  Dynnikov dyn(n);
  while (m--) {
    int i, s;
    scanf("%d%d", &i, &s);
    dyn.add(i, s);
  }
  dyn.print();
}

int main(int argc, char **argv) {
  if (argc > 1 && argv[1][0] == 'c') {
    convert();
  } else {
    test();
  }
  return 0;
}
