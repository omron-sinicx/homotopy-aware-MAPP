// A* search for commbined configuration space with homotopical information

#include <algorithm>
#include <vector>

#include "virtual_braid.hpp"
#include <cassert>
#include <cmath>

class AnonymousBhattacharya {
public:
  struct Point {
    int x, y;

    bool operator<(const Point &p2) const {
      if (this->x != p2.x) {
        return this->x < p2.x;
      }
      return this->y < p2.y;
    }

    bool operator>(const Point &p2) const { return p2 < *this; }
  };

  bool calc_ccw(const Point &p0, const Point &p1, const Point &p2) const {
    return (p1.x - p0.x) * (p2.y - p0.y) - (p2.x - p0.x) * (p1.y - p0.y) > 0;
  }

  // problem instance
  int noa, nov;
  std::vector<Point> coordinates;
  unsigned long long add_count = 0;

  // VirtualBraidPtr calc_next_braid(const VirtualBraidPtr &current_braid, const
  // int i,
  //                      const int target, std::vector<int> &next_pos) const;
  void inner_calc_next_braid(VirtualBraidPtr &current_braid, const int i,
                             const int target, std::vector<int> &next_pos);
};
