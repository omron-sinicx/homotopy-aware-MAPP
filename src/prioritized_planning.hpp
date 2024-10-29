#include "AnonymousBhattacharya.hpp"
#include <chrono>
#include <memory>

enum braid_types { generator_type, dynnikov_type };

class PrioritizedPlanning : public AnonymousBhattacharya {
public:
  braid_types braid_type = dynnikov_type;
  VirtualBraidPtr get_initial_braid(const int N = 0) const;
  VirtualBraidPtr get_braid_copy(const VirtualBraidPtr &b) const;

  struct Plan {
    int cost;
    std::vector<std::vector<int>> routes;
    int makespan;
    VirtualBraidPtr braid;
  };

  std::vector<std::vector<int>> edges, rev_edge_ids, edge_vertices;
  std::vector<int> starts, goals, obpts, edge_vertices_targets;
  int want;
  int count_closed_nodes = 0, count_all_nodes = 0;

  void reorder(std::vector<int> &order, std::vector<int> &rev_order,
               const int target, const int i,
               const std::vector<int> &pos) const {
    int rev = rev_order[i];
    for (int j = i - 1; j >= 0 && coordinates[target] < coordinates[pos[j]];
         j--) {
      int k = rev_order[j];
      order[k]++;
      order[rev]--;
      std::swap(rev_order[j], rev_order[j + 1]);
    }
    for (int j = i + 1;
         j < pos.size() && coordinates[target] > coordinates[pos[j]]; j++) {
      int k = rev_order[j];
      order[k]--;
      order[rev]++;
      std::swap(rev_order[j], rev_order[j - 1]);
    }
  }

  int heuristics(const int i, const int j) const {
    return (int)(abs(coordinates[i].x - coordinates[j].x) + 0.5) +
           (int)(abs(coordinates[i].y - coordinates[j].y) + 0.5);
  }

  std::vector<Plan> single_search(const std::vector<Plan> plans,
                                  const int agent_id);

  static int duration_to_int(const std::chrono::duration<double> &duration) {
    return std::chrono::duration_cast<std::chrono::milliseconds>(duration)
        .count();
  }

  std::vector<Plan> search(const std::string &log_filename = "",
                           const int time_limit = 0);
  void read_grid_map(int &n, int &m);
};
