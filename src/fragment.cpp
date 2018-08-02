#include "src/fragment.h"

Fragment::Fragment(const std::string& _chrom, const int32_t& _start, const size_t& _length) {
  chrom = _chrom;
  start = _start;
  length = _length;
}

Fragment::~Fragment() {}
