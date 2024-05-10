#ifndef VIRTUAL_BRAID_HEADER
#define VIRTUAL_BRAID_HEADER

#include <memory>
#include <tuple>
#include <vector>

class VirtualBraid {
public:
  int N;
  virtual void add(const int i, const bool sign) = 0;
  virtual bool operator<(const VirtualBraid &b) const = 0;
  virtual void post_process() = 0;
  virtual void print() const = 0;
};

using VirtualBraidPtr = std::shared_ptr<VirtualBraid>;

#endif
