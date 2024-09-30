#include <cstdio>
#include <cstring>
#include <queue>

const int dx[8] = {1, 1, 1, 0, 0, -1, -1, -1},
          dy[8] = {1, 0, -1, 1, -1, 1, 0, -1};

int main(int argc, char **argv) {
  int m, n;
  FILE *f = stdin;
  fscanf(f, "type octile\nheight %d\nwidth %d\nmap\n", &m, &n);
  std::vector<std::vector<bool>> map(n, std::vector<bool>(m));
  for (int i = 0; i < m; i++) {
    char s[1000];
    fscanf(f, "%s", s);
    for (int j = 0; j < n; j++) {
      map[j][m - i - 1] = (s[j] == '.' || s[j] == 'G');
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
  std::vector<std::vector<int>> e(nov);
  for (int x = 0; x < n; x++) {
    for (int y = 0; y < m; y++) {
      int v = m * x + y, id = ids[v];
      if (id == -1) {
        continue;
      }
      e[id].push_back(id);
      if (0 < x && ids[v - m] != -1) {
        e[id].push_back(ids[v - m]);
      }
      if (x + 1 < n && ids[v + m] != -1) {
        e[id].push_back(ids[v + m]);
      }
      if (0 < y && ids[v - 1] != -1) {
        e[id].push_back(ids[v - 1]);
      }
      if (y + 1 < m && ids[v + 1] != -1) {
        e[id].push_back(ids[v + 1]);
      }
    }
  }
  int obs_count = 0;
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
        obs_count++;
      }
    }
  }
  printf("%d\n", obs_count);
  return 0;
}
