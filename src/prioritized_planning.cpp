#include "prioritized_planning.hpp"
#include "braid.hpp"
#include "dynnikov.hpp"
#include <map>
#include <numeric>
#include <queue>
#include <set>

VirtualBraidPtr PrioritizedPlanning::get_initial_braid(const int N) const {
  return braid_type == generator_type ? VirtualBraidPtr(new Braid(N))
                                      : VirtualBraidPtr(new Dynnikov(N));
}

VirtualBraidPtr
PrioritizedPlanning::get_braid_copy(const VirtualBraidPtr &b) const {
  return braid_type == generator_type ? VirtualBraidPtr(new Braid(*b))
                                      : VirtualBraidPtr(new Dynnikov(*b));
}

std::vector<PrioritizedPlanning::Plan>
PrioritizedPlanning::single_search(const std::vector<Plan> plans,
                                   const int agent_id) {

  const int INF = 1000000000;

  struct Configure {
    int pos, time;
    VirtualBraidPtr braid;
    int base_plan_id;
    Configure() {}
    Configure(const int _pos, const int _time, const VirtualBraidPtr _braid,
              const int _base_plan_id)
        : pos(_pos), time(_time), braid(_braid), base_plan_id(_base_plan_id) {}
    bool operator<(const Configure &conf) const {
      if (pos != conf.pos) {
        return pos < conf.pos;
      } else if (time != conf.time) {
        return time < conf.time;
      } else if (base_plan_id != conf.base_plan_id) {
        return base_plan_id < conf.base_plan_id;
      }
      return *braid < *conf.braid;
    }
  };
  struct Node {
    Configure conf;
    std::vector<int> order, rev_order, pos;
    int dist;
    int back, back_edge;
    bool open;
  };

  int start = starts[agent_id];
  int goal = goals[agent_id];

  int noo = obpts.size();
  std::vector<int> combined_starts = obpts;
  combined_starts.insert(combined_starts.end(), starts.begin(), starts.end());
  std::vector<int> start_order(noo + agent_id + 1),
      start_rev_order(noo + agent_id + 1), start_pos(noo + agent_id + 1);
  std::iota(start_rev_order.begin(), start_rev_order.end(), 0);
  std::vector<Point> start_coordinates(noo + agent_id + 1);
  for (int i = 0; i <= noo + agent_id; i++) {
    start_coordinates[i] = coordinates[combined_starts[i]];
  }
  std::sort(start_rev_order.begin(), start_rev_order.end(),
            [&start_coordinates](const int i, const int j) {
              return start_coordinates[i] < start_coordinates[j];
            });
  for (int i = 0; i <= noo + agent_id; i++) {
    start_order[start_rev_order[i]] = i;
  }
  for (int i = 0; i <= noo + agent_id; i++) {
    start_pos[i] = combined_starts[start_rev_order[i]];
  }

  std::vector<Node> nodes;
  using QPair = std::pair<std::pair<int, int>, int>;
  std::priority_queue<QPair, std::vector<QPair>, std::greater<QPair>> Q;
  std::map<Configure, int> node_index;
  std::vector<int> goal_time(plans.size(), 0);

  std::vector<std::vector<std::vector<int>>> dist_to_go(plans.size());

  std::vector<std::vector<std::vector<bool>>> plan_cells(plans.size());
  std::vector<std::vector<std::vector<std::vector<bool>>>> plan_edge_cells(plans.size());
  
  for (int base_plan_id = 0; base_plan_id < plans.size(); base_plan_id++) {
    auto &plan = plans[base_plan_id];

    int height = plan.makespan + 1;
    std::vector<std::vector<bool>> cells(edges.size(),
                                         std::vector<bool>(height, true));
    std::vector<std::vector<std::vector<bool>>> edge_cells(edges.size(),
							   std::vector<std::vector<bool>>(height));
    for(int v = 0; v < edges.size(); v++){
      for(int t = 0; t < height; t++){
	edge_cells[v][t] = std::vector<bool>(edges[v].size(), true);
      }
    }
    std::vector<int> cur = starts;
    cur.resize(agent_id);
    for (int t = 0; t < height; t++) {
      for (int i = 0; i < agent_id; i++) {
        cells[cur[i]][t] = false;
        if (cur[i] == goal) {
          goal_time[base_plan_id] = t + 1;
        }
        if (plan.routes[i].size() <= t) {
          continue;
        }
        int next_cur = edges[cur[i]][plan.routes[i][t]];
	int rev_edge_id = rev_edge_ids[cur[i]][plan.routes[i][t]];
	edge_cells[next_cur][t][rev_edge_id] = false;
        cur[i] = next_cur;
      }
    }

    for (int i = agent_id + 1; i < starts.size(); i++) {
      if (starts[i] == goal) {
	continue;
      }
      for (int t = 0; t < height; t++) {
        cells[starts[i]][t] = false;
      }
    }

    plan_cells[base_plan_id] = cells;
    plan_edge_cells[base_plan_id] = edge_cells;

    auto &dist = dist_to_go[base_plan_id];
    dist = std::vector<std::vector<int>>(edges.size(),
                                         std::vector<int>(height, INF));
    std::queue<std::pair<int, int>> bfs;
    for (int t = goal_time[base_plan_id]; t < height; t++) {
      bfs.push(std::make_pair(goal, t));
      dist[goal][t] = 0;
    }
    while (!bfs.empty()) {
      int i = bfs.front().first, t = bfs.front().second;
      bfs.pop();
      if (t == height - 1) {
        for (auto j : edges[i]) {
          if (cells[j][t] && dist[j][t] == INF) {
            dist[j][t] = dist[i][t] + 1;
            bfs.push(std::make_pair(j, t));
          }
        }
      }
      if (t > 0) {
        for (int e = 0 ; e < edges[i].size(); e++) {
	  int j = edges[i][e];
          if (cells[j][t - 1] && edge_cells[j][t-1][rev_edge_ids[i][e]] && dist[j][t - 1] == INF) {
            dist[j][t - 1] = dist[i][t] + 1;
            bfs.push(std::make_pair(j, t - 1));
          }
        }
      }
    }

    if (dist[start][0] == INF) {
      continue;
    }

    Node node;
    node.conf = Configure(start, 0, get_initial_braid(noo + agent_id + 1),
                          base_plan_id);
    node.order = start_order;
    node.rev_order = start_rev_order;
    node.pos = start_pos;
    node.dist = plan.cost;
    node.back = node.back_edge = -1;
    node.open = true;

    int id = nodes.size();
    nodes.push_back(node);
    int h = dist[start][0];
    Q.push(std::make_pair(std::make_pair(node.dist + h, h), id));
    node_index[node.conf] = id;
  }

  if (nodes.size() == 0) {
    char message[100];
    sprintf(message, "Path Planning Failed at Agent %d\n", agent_id);
    throw std::runtime_error(std::string(message));
  }
  std::vector<Plan> new_plans;

  int number_of_plans = want;

  int max_height = 0;
  for (auto &plan : plans) {
    max_height = std::max(max_height, plan.makespan + 1);
  }
  // std::vector<std::vector<int>> braid_count(coordinates.size(),
  // std::vector<int>(max_height,0));

  count_closed_nodes = 0;

  while (!Q.empty() && new_plans.size() < number_of_plans) {
    QPair qpair = Q.top();
    Q.pop();
    int id = qpair.second;
    if (!nodes[id].open)
      continue;
    nodes[id].open = false;
    count_closed_nodes++;

    Node node = nodes[id];
    int node_dist = node.dist;
    Configure &conf = node.conf;
    int t = conf.time;

    // if(braid_count[conf.pos][t] >= number_of_plans){
    // continue;
    //}
    // braid_count[conf.pos][t]++;

    auto &plan = plans[conf.base_plan_id];

    if (conf.pos == goal && t == plan.makespan) {
      Plan new_plan;
      std::vector<int> route;
      int cur = id;
      while (nodes[cur].back != -1 && nodes[nodes[cur].back].conf.pos == goal) {
        cur = nodes[cur].back;
      }
      while (nodes[cur].back != -1) {
        route.push_back(nodes[cur].back_edge);
        cur = nodes[cur].back;
      }
      std::reverse(route.begin(), route.end());
      new_plan.routes = plan.routes;
      new_plan.routes.push_back(route);
      new_plan.cost = node_dist;
      new_plan.makespan = std::max(plan.makespan, (int)route.size());
      new_plan.braid = conf.braid;
      new_plans.push_back(new_plan);
      continue;
    }
    // move determined agents according to plan
    std::vector<int> order = node.order, rev_order = node.rev_order,
                     pos = node.pos;
    VirtualBraidPtr braid = get_braid_copy(conf.braid);
    int s = pos[order[noo + agent_id]];
    for (int i = 0; i < agent_id; i++) {
      if (plan.routes[i].size() <= t) {
        continue;
      }
      int order_i = order[noo + i], source = pos[order_i];
      int target = edge_vertices[source][plan.routes[i][t]];
      reorder(order, rev_order, target, order_i, pos);
      inner_calc_next_braid(braid, order_i, target, pos);
    }
    braid->post_process();
    for (int edge_id = 0; edge_id < edges[s].size(); edge_id++) {
      int target = edges[s][edge_id], target_time = std::min(t + 1, plan.makespan);
      if(!plan_cells[conf.base_plan_id][target][target_time] || !plan_edge_cells[conf.base_plan_id][s][t][edge_id]){
	continue;
      }
      Configure new_conf;
      new_conf.pos = target;
      new_conf.time = target_time;
      new_conf.base_plan_id = conf.base_plan_id;
      
      int h = dist_to_go[new_conf.base_plan_id][new_conf.pos][new_conf.time];
      if (h == INF) {
        continue;
      }

      // calculate new braid
      new_conf.braid = get_braid_copy(braid);
      auto new_order = order;
      auto new_rev_order = rev_order;
      auto new_pos = pos;
      int edge_vertex = edge_vertices[s][edge_id], order_i = new_order[noo + agent_id];
      reorder(new_order, new_rev_order, edge_vertex, order_i, new_pos);
      inner_calc_next_braid(new_conf.braid, order_i, edge_vertex, new_pos);
      for (int i = 0; i <= agent_id; i++) {
	if (i < agent_id && plan.routes[i].size() <= t) {
	  continue;
	}
	int order_i = new_order[noo + i];
	int target = edge_vertices_targets[new_pos[order_i]];
	assert(target!=-1);
	reorder(new_order, new_rev_order, target, order_i, new_pos);
	inner_calc_next_braid(new_conf.braid, order_i, target, new_pos);
      }
      new_conf.braid->post_process();

      int new_dist = node_dist + (s == goal && target == goal &&
                                          t >= goal_time[conf.base_plan_id]
                                      ? 0
                                      : 1);
      if (node_index.find(new_conf) == node_index.end()) {
        Node new_node;
        new_node.conf = new_conf;
        new_node.dist = new_dist;
        new_node.back = id;
        new_node.back_edge = edge_id;
        new_node.open = true;
        new_node.order = new_order;
        new_node.rev_order = new_rev_order;
        new_node.pos = new_pos;
        int new_id = nodes.size();
        node_index[new_conf] = new_id;
        nodes.push_back(new_node);
        Q.push(std::make_pair(std::make_pair(new_dist + h, h), new_id));
      } else {
        int old_id = node_index[new_conf];
        auto &old_node = nodes[old_id];
        if (old_node.dist > new_dist) {
          old_node.dist = new_dist;
          old_node.back = id;
          old_node.back_edge = edge_id;
          Q.push(std::make_pair(std::make_pair(new_dist + h, h), old_id));
        }
      }
    }
  }
  count_all_nodes = nodes.size();
  return new_plans;
}

std::vector<PrioritizedPlanning::Plan>
PrioritizedPlanning::search(const std::string &log_filename,
                            const int time_limit) {
  std::vector<Plan> plans(1);
  plans[0].routes = std::vector<std::vector<int>>(0);
  plans[0].cost = 0;
  plans[0].makespan = 0;
  plans[0].braid = get_initial_braid(obpts.size());

  FILE *log_file = NULL;
  if (log_filename != "") {
    if (log_filename == std::string("stderr")) {
      log_file = stderr;
    } else {
      log_file = fopen(log_filename.c_str(), "w");
    }
    assert(log_file != NULL);
  }
  auto start_time = std::chrono::system_clock::now();
  for (int i = 0; i < noa; i++) {
    plans = single_search(plans, i);
    if (log_file != NULL) {
      int runtime =
	duration_to_int(std::chrono::system_clock::now() - start_time);
      fprintf(log_file, "%d %d %d %llu %d %d ", runtime, plans[0].makespan,
              plans[plans.size() - 1].makespan, add_count, count_closed_nodes,
              count_all_nodes);
      if (braid_type == generator_type) {
        auto b = dynamic_cast<Braid &>(*plans[0].braid);
        auto bl = dynamic_cast<Braid &>(*plans[plans.size() - 1].braid);
        fprintf(log_file, "%d %d\n", (int)b.word.size(), (int)bl.word.size());
      } else {
        auto b = dynamic_cast<Dynnikov &>(*plans[0].braid);
        auto bl = dynamic_cast<Dynnikov &>(*plans[plans.size() - 1].braid);
        gmp_fprintf(log_file, "%Zd %Zd\n", b.get_max_abs_cd(),
                    bl.get_max_abs_cd());
      }
      if (time_limit > 0 && runtime > time_limit) {
        break;
      }
      fflush(log_file);
    }
  }
  if (log_file != NULL) {
    fclose(log_file);
  }
  return plans;
}

const int dx[8] = {1, 1, 1, 0, 0, -1, -1, -1},
          dy[8] = {1, 0, -1, 1, -1, 1, 0, -1};

void PrioritizedPlanning::read_grid_map(int &n, int &m) {
  std::vector<std::vector<bool>> map;
  char map_type[100];
  scanf("%s", map_type);
  if (!strcmp(map_type, "empty")) {
    scanf("%d%d", &n, &m);
    map = std::vector<std::vector<bool>>(n, std::vector<bool>(m, true));
  } else {
    FILE *f = fopen(map_type, "r");
    if (f == NULL) {
      fprintf(stderr, "Invalid map file");
      exit(1);
    }
    fscanf(f, "type octile\nheight %d\nwidth %d\nmap\n", &m, &n);
    map = std::vector<std::vector<bool>>(n, std::vector<bool>(m));
    assert(n<=1000);
    for (int i = 0; i < m; i++) {
      char s[1000];
      fscanf(f, "%s", s);
      for (int j = 0; j < n; j++) {
        map[j][m - i - 1] = (s[j] == '.' || s[j] == 'G');
      }
    }
  }
  int nov = 0;
  std::vector<int> ids(n * m);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      if (map[i][j]) {
        ids[m * i + j] = nov++;
      } else {
        ids[m * i + j] = -1;
      }
    }
  }
  edges = std::vector<std::vector<int>>(nov);
  for (int x = 0; x < n; x++) {
    for (int y = 0; y < m; y++) {
      int v = m * x + y, id = ids[v];
      if (id == -1) {
        continue;
      }
      edges[id].push_back(id);
      if (0 < x && ids[v - m] != -1) {
        edges[id].push_back(ids[v - m]);
      }
      if (x + 1 < n && ids[v + m] != -1) {
        edges[id].push_back(ids[v + m]);
      }
      if (0 < y && ids[v - 1] != -1) {
        edges[id].push_back(ids[v - 1]);
      }
      if (y + 1 < m && ids[v + 1] != -1) {
        edges[id].push_back(ids[v + 1]);
      }
    }
  }
  coordinates = std::vector<AnonymousBhattacharya::Point>(nov);
  for (int x = 0; x < n; x++) {
    for (int y = 0; y < m; y++) {
      int id = ids[m * x + y];
      if (id != -1) {
        coordinates[id].x = 2*x;
        coordinates[id].y = 2*y;
      }
    }
  }
  obpts = std::vector<int>();
  std::vector<std::vector<bool>> visit(n, std::vector<bool>(m, false));
  for (int x = 0; x < n; x++) {
    for (int y = 0; y < m; y++) {
      if (map[x][y] or visit[x][y]) {
        continue;
      }
      std::queue<std::pair<int, int>> Q;
      Q.push(std::make_pair(x, y));
      visit[x][y] = true;
      bool outer = false;
      while (!Q.empty()) {
        auto p = Q.front();
        Q.pop();
        if (p.first == 0 || p.first == n - 1 || p.second == 0 ||
            p.second == m - 1) {
          outer = true;
        }
        for (int v = 0; v < 8; v++) {
          int nx = p.first + dx[v], ny = p.second + dy[v];
          if (0 <= nx && nx < n && 0 <= ny && ny < m && !map[nx][ny] &&
              !visit[nx][ny]) {
            Q.push(std::make_pair(nx, ny));
            visit[nx][ny] = true;
          }
        }
      }
      if (!outer) {
        obpts.push_back(coordinates.size());
        coordinates.push_back(AnonymousBhattacharya::Point(2*x, 2*y));
      }
    }
  }
  int k;
  scanf("%d", &k);
  std::vector<int> x0(k), y0(k), x1(k), y1(k);
  for (int i = 0; i < k; i++) {
    scanf("%d%d%d%d", &x0[i], &y0[i], &x1[i], &y1[i]);
  }
  noa = k;
  starts = std::vector<int>(k);
  goals = std::vector<int>(k);
  for (int i = 0; i < k; i++) {
    starts[i] = ids[m * x0[i] + y0[i]];
    if (starts[i] == -1) {
      fprintf(stderr, "Invalid start for agent %d\n", i);
      exit(1);
    }
    goals[i] = ids[m * x1[i] + y1[i]];
    if (goals[i] == -1) {
      fprintf(stderr, "Invalid goal for agent %d\n", i);
      exit(1);
    }
  }

  rev_edge_ids = std::vector<std::vector<int>>(edges.size());
  edge_vertices = std::vector<std::vector<int>>(edges.size());
  edge_vertices_targets = std::vector<int>(coordinates.size(), -1);
  for(int i = 0; i < edges.size(); i++){
    rev_edge_ids[i] = std::vector<int>(edges[i].size());
    edge_vertices[i] = std::vector<int>(edges[i].size());
    for(int e = 0 ; e < edges[i].size(); e++){
      int j = edges[i][e];
      for(int ej =0 ; ej < edges[j].size(); ej++){
	if(edges[j][ej] == i){
	  rev_edge_ids[i][e] = ej;
	}
      }
      edge_vertices[i][e] = coordinates.size();
      edge_vertices_targets.push_back(j);
      coordinates.push_back(Point((coordinates[i].x+coordinates[j].x)/2,
				  (coordinates[i].y+coordinates[j].y)/2));
    }
  }
}
