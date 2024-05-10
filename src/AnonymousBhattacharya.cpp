#include "AnonymousBhattacharya.hpp"

void AnonymousBhattacharya::inner_calc_next_braid(
    VirtualBraidPtr &current_braid, const int i, const int target,
    std::vector<int> &next_pos) {
  // return the braid after agent i move to target

  int start = next_pos[i];
  bool dir = (coordinates[start] < coordinates[target]);
  next_pos[i] = target;
  for (int j = dir ? i + 1 : i - 1; dir ? j < next_pos.size() : j >= 0;
       dir ? j++ : j--) {
    int pos_j = next_pos[j];
    if (dir ? !(coordinates[target] > coordinates[pos_j])
            : !(coordinates[pos_j] > coordinates[target]))
      break;
    std::swap(next_pos[dir ? j - 1 : j + 1], next_pos[j]);
    bool ccw =
        calc_ccw(coordinates[start], coordinates[pos_j], coordinates[target]);
    current_braid->add(dir ? j - 1 : j, ccw);
    add_count++;
  }
}

/*VirtualBraidPtr AnonymousBhattacharya::calc_next_braid(const VirtualBraidPtr
  &current_braid, const int i, const int target, std::vector<int> &next_pos)
  const { VirtualBraidPtr new_braid = std::make_unique(*current_braid);
  inner_calc_next_braid(new_braid, i, target, next_pos);
  new_braid.post_process();
  return new_braid;
  return current_braid;
  }*/
