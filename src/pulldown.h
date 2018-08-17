#ifndef SRC_PULLDOWN_H__
#define SRC_PULLDOWN_H__

#include "src/bingenerator.h"
#include "src/fragment.h"
#include "src/options.h"
#include "src/peak_intervals.h"

#include <vector>


class Pulldown {
 public:
  Pulldown(const Options& options, const GenomeBin& gbin);
  void Perform(vector<Fragment>* output_fragments, PeakIntervals* pintervals);

 private:
  std::string chrom;
  int32_t start;
  int32_t end;
  int numcopies;
  float gamma_alpha, gamma_beta;
  float ratio_beta;
  bool debug_pulldown;
};

#endif  // SRC_PULLDOWN_H__
