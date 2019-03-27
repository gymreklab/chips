#ifndef SRC_PEAKINTERVALS_H__
#define SRC_PEAKINTERVALS_H__

#include <map>
#include <vector>
#include "src/fragment.h"
#include "src/peak_loader.h"
#include "src/options.h"

class PeakIntervals {
  /*
    This class contains all peaks and their scores
    It is used to determine whether a fragment overlaps a peak and what the score is
   */
 public:
  PeakIntervals(const Options& options, const std::string peakfile, const std::string peakfileType, const std::string bamfile,
		const std::int32_t count_colidx);
  virtual ~PeakIntervals();

  /* Get score of peak overlapping fragment */
  float GetOverlap(const Fragment& frag, int& peakIndexStart);
  void resetSearchScope(const int index);
  float total_bound_length;
  float total_genome_length;

 private:
  // peakmap:  key: chromID,  data: fragments
  std::map<std::string, std::vector<Fragment> > peak_map;
  /* Load peaks from file */
  bool LoadPeaks(const Options& options, const std::string peakfile, const std::string peakfileType, const std::string bamfile,
		 const std::int32_t count_colidx);
  float SearchList(const Fragment& frag, int& peakIndexStart);
};

#endif  // SRC_PEAKINTERVALS_H__
