#include "AnonymousBhattacharya.hpp"
#include <chrono>

class PrioritizedPlanning : public AnonymousBhattacharya {
public:
  struct Plan {
    int cost;
    std::vector<std::vector<int>> routes;
    int makespan;
    Braid braid;
  };

  std::vector<std::vector<int>> edges;
  std::vector<int> starts, goals;
  int want;

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
                                  const int agent_id) const;

  static int duration_to_int(const std::chrono::duration<double> &duration) {
    return std::chrono::duration_cast<std::chrono::milliseconds>(duration)
        .count();
  }

  std::vector<Plan> search(const std::string &log_filename = "") const {
    std::vector<Plan> plans(1);
    plans[0].routes = std::vector<std::vector<int>>(0);
    plans[0].cost = 0;
    plans[0].makespan = 0;
    plans[0].braid = Braid(0);

    FILE *log_file = NULL;
    if(log_filename != ""){
      log_file = fopen(log_filename.c_str(), "w");
    }
    auto start_time = std::chrono::system_clock::now();
    for (int i = 0; i < noa; i++) {
      plans = single_search(plans, i);
      if(log_file != NULL){
	int runtime = duration_to_int(std::chrono::system_clock::now() - start_time);
	fprintf(log_file, "%d %d %d\n", runtime, (int)plans[0].braid.word.size(), (int)plans[0].makespan);
      }
    }
    if (log_file != NULL) {
      fclose(log_file);
    }
    return plans;
  }
};
