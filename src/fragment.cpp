#include "src/fragment.h"

Fragment::Fragment(const std::string& _chrom, const int32_t& _start, const size_t& _length, const float& _score){
  chrom = _chrom;
  start = _start;
  length = _length;
  score = _score;
  orig_score = _score;
}

Fragment::~Fragment() {}
