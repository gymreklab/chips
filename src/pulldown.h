#ifndef SRC_PULLDOWN_H__
#define SRC_PULLDOWN_H__

#include "src/bingenerator.h"
#include "src/fragment.h"
#include "src/options.h"
#include "src/peak_intervals.h"

#include <vector>


class Pulldown {
 public:
  Pulldown(const Options& options, const GenomeBin& gbin,
            std::string& _prev_chrom, int& _peakIndexStart, int& _start_offset);
  void Perform(vector<Fragment>* output_fragments, PeakIntervals* pintervals);

 private:
  std::string chrom;
  std::int32_t start;
  std::int32_t end;
  int numcopies;
  float gamma_alpha, gamma_beta;
  float ratio_beta;
  bool debug_pulldown;

  int peakIndexStart;
  std::string prev_chrom;
  int start_offset;
};
#endif  // SRC_PULLDOWN_H__
