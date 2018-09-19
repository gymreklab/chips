#ifndef SRC_PEAKINTERVALS_H__
#define SRC_PEAKINTERVALS_H__

#include <map>
#include <vector>
#include "src/fragment.h"
#include "src/peak_loader.h"

class PeakIntervals {
  /*
    This class contains all peaks and their scores
    It is used to determine whether a fragment overlaps a peak and what the score is
   */
 public:
  PeakIntervals(const std::string peakfile, const std::string peakfileType, const std::string bamfile, const std::int32_t count_colidx);
  virtual ~PeakIntervals();

  /* Get score of peak overlapping fragment */
  float GetOverlap(const Fragment& frag);
  void resetSearchScope(const int index);
  int peakIndexStart;
  double prob_pd_given_b;

 private:
  // peakmap:  key: chromID,  data: fragments
  std::map<std::string, std::vector<Fragment> > peak_map;
  std::uint32_t max_coverage;
  /* Load peaks from file */
  bool LoadPeaks(const std::string peakfile, const std::string peakfileType, const std::string bamfile, const std::int32_t count_colidx);
  float SearchList(const Fragment& frag);
};

#endif  // SRC_PEAKINTERVALS_H__
