#include <tuple>
#include <vector>

class Braid {
public:
  int N;
  bool sign;
  int main_i;
  Braid() {}
  Braid(const int _N) : N(_N) {
    sign = false;
    main_i = N;
  }
  std::vector<std::pair<int, bool>> word;
  void reduce_freely();
  void reduce_fully();
  void add(const int i, const bool sign);
  std::pair<int, bool> get_main() const;
  bool is_positive() const;

  bool operator<(const Braid &b) const {
    if (sign != b.sign)
      return sign < b.sign;
    if (main_i != b.main_i)
      return sign ^ (main_i < b.main_i);
    Braid diff;
    diff.N = N;
    diff.word.resize(word.size() + b.word.size());
    for (int i = 0; i < word.size(); i++) {
      auto w = word[word.size() - i - 1];
      diff.word[i] = std::make_pair(w.first, !w.second);
    }
    std::copy(b.word.begin(), b.word.end(), diff.word.begin() + word.size());
    diff.reduce_fully();
    return diff.is_positive();
  }

  void print() const {
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
