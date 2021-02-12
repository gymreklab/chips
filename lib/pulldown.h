#ifndef SRC_PULLDOWN_H__
#define SRC_PULLDOWN_H__

#include "bingenerator.h"
#include "fragment.h"
#include "options.h"
#include "peak_intervals.h"

#include <vector>
#include <random>


class Pulldown {
 public:
  Pulldown(const Options& options, const GenomeBin& gbin,
            std::string& _prev_chrom, int& _peakIndexStart, int& _start_offset);
  void Perform(vector<Fragment>* output_fragments, PeakIntervals* pintervals, std::mt19937& rng);

 private:
  std::string chrom;
  std::int32_t start;
  std::int32_t end;
  int numcopies;
  float gamma_k, gamma_theta;
  float ratio_beta;
  bool debug_pulldown;

  int peakIndexStart;
  std::string prev_chrom;
  int* start_offset_ptr;
  unsigned seed;
};
#endif  // SRC_PULLDOWN_H__
