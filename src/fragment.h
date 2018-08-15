#ifndef SRC_FRAGMENT_H__
#define SRC_FRAGMENT_H__

#include <stdint.h>
#include <string>

class Fragment {
 public:
  Fragment(const std::string& _chrom, const int32_t& _start, const size_t& _length, const float& _score=0);
  virtual ~Fragment();
  
  std::string chrom;
  int32_t start;
  size_t length;
  float score; // store tagcounts or probabilty of being bound
};

#endif  // SRC_FRAGMENT_H__
