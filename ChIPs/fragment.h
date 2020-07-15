#ifndef SRC_FRAGMENT_H__
#define SRC_FRAGMENT_H__

#include <stdint.h>
#include <string>

class Fragment {
 public:
  Fragment(const std::string& _chrom, const std::int32_t& _start, const std::size_t& _length, const float& _score=0);
  virtual ~Fragment();
  
  std::string chrom;
  std::int32_t start;
  std::size_t length;
  float score; // store tagcounts or probabilty of being bound
  float orig_score; // store score before scaling
};

#endif  // SRC_FRAGMENT_H__
