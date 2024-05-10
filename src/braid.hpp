#include "virtual_braid.hpp"
#include <tuple>
#include <vector>

class Braid : public VirtualBraid {
public:
  Braid() {}
  Braid(const int _N) {
    N = _N;
    sign = false;
    main_i = N;
  }
  Braid(const VirtualBraid &_b);
  void add(const int i, const bool sign) override;
  bool operator<(const VirtualBraid &b) const override;
  void post_process() override { reduce_fully(); }
  bool sign;
  int main_i;
  std::vector<std::pair<int, bool>> word;
  void reduce_freely();
  void reduce_fully();
  std::pair<int, bool> get_main() const;
  bool is_positive() const;

  void print() const override {
    printf("%d ", (int)word.size());
    for (int i = 0; i < word.size(); i++) {
      printf("%d %c ", word[i].first, word[i].second ? '+' : '-');
    }
    putchar('\n');
  }
};

class PureBraid {
public:
  int N;
  std::vector<std::tuple<int, int, bool>> word;
  PureBraid() {}
  PureBraid(const int _N) : N(_N) {}
  void add(const int i, const int j, const bool sign);

  Braid reduced, rev_reduced;
  bool sign;
  int main_i;

  void calc_reduced();

  bool operator==(const PureBraid &b) const {
    if (sign != b.sign)
      return false;
    if (main_i != b.main_i)
      return false;
    Braid diff = rev_reduced;
    diff.word.insert(diff.word.end(), b.reduced.word.begin(),
                     b.reduced.word.end());
    diff.reduce_fully();
    return diff.word.size() == 0;
  }

  bool operator<(const PureBraid &b) const {
    if (sign != b.sign)
      return sign < b.sign;
    if (main_i != b.main_i)
      return sign ^ (main_i < b.main_i);
    Braid diff = rev_reduced;
    diff.word.insert(diff.word.end(), b.reduced.word.begin(),
                     b.reduced.word.end());
    diff.reduce_fully();
    return diff.is_positive();
  }

  void print() {
    printf("%d ", (int)word.size());
    for (auto &t : word) {
      printf("(%d %d)%c ", std::get<0>(t), std::get<1>(t),
             std::get<2>(t) ? '+' : '-');
    }
    putchar('\n');
  }
};
