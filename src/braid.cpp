#include "braid.hpp"
#include <algorithm>
#include <cassert>

void Braid::reduce_freely() {
  int resize = 0;
  for (int i = 0; i < word.size(); i++) {
    if (resize > 0 && word[resize - 1].first == word[i].first &&
        word[resize - 1].second != word[i].second) {
      resize--;
    } else {
      word[resize++] = word[i];
    }
  }
  word.resize(resize);
}

void Braid::reduce_fully() {
  while (true) {
    reduce_freely();
    std::vector<int> memo(N, -1);
    int a = -1, b = -1;
    for (int i = 0; i < word.size(); i++) {
      int j = word[i].first;
      if (memo[j] != -1 && word[memo[j]].second != word[i].second &&
          (j == 0 || memo[j - 1] < memo[j])) {
        a = memo[j];
        b = i;
        break;
      }
      memo[j] = i;
    }
    if (a == -1) {
      if (word.size() == 0) {
        sign = false;
        main_i = N;
      } else {
        for (int i = 0; i < N; i++) {
          if (memo[i] != -1) {
            sign = word[memo[i]].second;
            main_i = i;
            break;
          }
        }
      }
      break;
    }
    int bundle_j = word[a].first;
    bool bundle_sign = word[a].second;
    std::vector<std::pair<int, bool>> new_word(a);
    std::copy(word.begin(), word.begin() + a, new_word.begin());
    for (int i = a + 1; i < b; i++) {
      if (word[i].first == bundle_j + 1) {
        new_word.push_back(std::make_pair(bundle_j + 1, !bundle_sign));
        new_word.push_back(std::make_pair(bundle_j, word[i].second));
        new_word.push_back(std::make_pair(bundle_j + 1, bundle_sign));
      } else {
        new_word.push_back(word[i]);
      }
    }
    new_word.insert(new_word.end(), word.begin() + b + 1, word.end());
    word = new_word;
  }
}

std::pair<int, bool> Braid::get_main() const {
  // assume this is fully reduced
  return std::make_pair(main_i, sign);
}

bool Braid::is_positive() const { return get_main().second; }

void Braid::add(const int i, const bool sign) {
  word.push_back(std::make_pair(i, sign));
}

void PureBraid::add(const int i, const int j, const bool sign) {
  assert(i < j);
  int last = word.size() - 1;
  if (last >= 0 && std::get<0>(word[last]) == i &&
      std::get<1>(word[last]) == j && std::get<2>(word[last]) != sign) {
    word.pop_back();
  } else {
    word.push_back(std::make_tuple(i, j, sign));
  }
}

void PureBraid::calc_reduced() {
  reduced.N = N;
  reduced.word.clear();
  for (auto &tpl : word) {
    int i = std::get<0>(tpl), j = std::get<1>(tpl);
    bool sign = std::get<2>(tpl);
    for (int k = j - 1; k > i; k--) {
      reduced.word.push_back(std::make_pair(k, true));
    }
    reduced.word.push_back(std::make_pair(i, sign));
    reduced.word.push_back(std::make_pair(i, sign));
    for (int k = i + 1; k < j; k++) {
      reduced.word.push_back(std::make_pair(k, false));
    }
  }
  reduced.reduce_fully();
  auto main_gen = reduced.get_main();
  main_i = main_gen.first;
  sign = main_gen.second;
  rev_reduced.N = N;
  rev_reduced.word.resize(reduced.word.size());
  for (int i = 0; i < reduced.word.size(); i++) {
    auto &rev = reduced.word[reduced.word.size() - i - 1];
    rev_reduced.word[i] = std::make_pair(rev.first, !rev.second);
  }
}
