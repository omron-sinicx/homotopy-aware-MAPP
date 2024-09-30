#include "braid.hpp"
#include "dynnikov.hpp"
#include <algorithm>
#include <cassert>

mpz_class get_max_cd(const Dynnikov &b) {
  mpz_class max_cd(0);
  for (auto &v : b.cd) {
    max_cd = std::max(max_cd, v > 0 ? v : -v);
  }
  return max_cd;
}

int main(int argc, char **argv) {
  int n = std::stoi(argv[1]), L = std::stoi(argv[2]), T = std::stoi(argv[3]),
      I = std::stoi(argv[4]);
  for (int t = 0; t < T; t++) {
    Dynnikov b(n);
    for (int l = 1; l <= L * I; l++) {
      int i = rand() % (n - 1);
      bool sign = rand() % 2;
      b.add(i, sign);
      if (l % L == 0) {
        auto max_cd = get_max_cd(b);
        gmp_printf("%.Zd\n", max_cd);
      }
    }
  }
  return 0;
}
