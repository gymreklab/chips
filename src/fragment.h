#ifndef SRC_FRAGMENT_H__
#define SRC_FRAGMENT_H__

#include <stdint.h>
#include <string>

class Fragment {
 public:
  Fragment(const std::string& _chrom, const int32_t& _start, const size_t& _length);
  virtual ~Fragment();
  
  std::string chrom;
  int32_t start;
  size_t length;
};

#endif  // SRC_FRAGMENT_H__
