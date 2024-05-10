#include "virtual_braid.hpp"
#include <gmpxx.h>
#include <iostream>
#include <vector>

mpz_class pls(const mpz_class &x) { return std::max(x, mpz_class(0)); }

mpz_class mns(const mpz_class &x) { return std::min(x, mpz_class(0)); }

class Dynnikov : public VirtualBraid {
public:
  std::vector<mpz_class> cd;

  Dynnikov() {}
  Dynnikov(const int _N) {
    N = _N;
    cd = std::vector<mpz_class>(2 * std::max(N, 1) - 2);
    for (int i = 0; i < cd.size(); i++) {
      cd[i] = -(i % 2);
    }
  }
  Dynnikov(const VirtualBraid &_b) {
    auto b = dynamic_cast<const Dynnikov &>(_b);
    N = b.N;
    cd = b.cd;
  }

  void add(const int i, const bool sign) override {
    if (i == 0) {
      if (sign) {
        mpz_class a = cd[0], b = cd[1];
        cd[0] = -b + pls(a + pls(b));
        cd[1] = a + pls(b);
      } else {
        mpz_class a = cd[0], b = cd[1];
        cd[0] = b - pls(pls(b) - a);
        cd[1] = pls(b) - a;
      }
    } else {
      mpz_class a = cd[2 * i - 2], b = cd[2 * i - 1], c = cd[2 * i],
                d = cd[2 * i + 1];
      if (sign) {
        mpz_class e = a + mns(b) - c - pls(d);
        cd[2 * i - 2] = a - pls(b) - pls(pls(d) + e);
        cd[2 * i - 1] = d + mns(e);
        cd[2 * i] = c - mns(d) - mns(mns(b) - e);
        cd[2 * i + 1] = b - mns(e);
      } else {
        mpz_class e = a - mns(b) - c + pls(d);
        cd[2 * i - 2] = a + pls(b) + pls(pls(d) - e);
        cd[2 * i - 1] = d - pls(e);
        cd[2 * i] = c + mns(d) + mns(mns(b) + e);
        cd[2 * i + 1] = b + pls(e);
      }
    }
  }

  bool operator<(const VirtualBraid &_b) const override {
    auto b = dynamic_cast<const Dynnikov &>(_b);
    for (int i = 0; i < cd.size(); i++) {
      if (cd[i] != b.cd[i])
        return cd[i] < b.cd[i];
    }
    return false;
  }

  bool operator==(const Dynnikov &b) const {
    for (int i = 0; i < cd.size(); i++) {
      if (cd[i] != b.cd[i])
        return false;
    }
    return true;
  }

  void post_process() override {}

  void print() const override {
    for (int i = 0; i < cd.size(); i++) {
      std::cout << cd[i] << ' ';
    }
    std::cout << '\n';
  }

  mpz_class get_max_abs_cd() const {
    mpz_class m(0);
    for (auto &v : cd) {
      m = std::max(m, v > 0 ? v : -v);
    }
    return m;
  }
};
