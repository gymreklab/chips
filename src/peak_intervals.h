#ifndef SRC_PEAKINTERVALS_H__
#define SRC_PEAKINTERVALS_H__

#include "src/fragment.h"

using namespace std;

class PeakIntervals {
  /*
    This class contains all peaks and their scores
    It is used to determine whether a fragment overlaps a peak and what the score is
   */
 public:
  PeakIntervals(const std::string& peakfile);
  virtual ~PeakIntervals();

  /* Get score of peak overlapping fragment */
  float GetOverlap(const Fragment& frag);

 private:
  /* Load peaks from file */
  bool LoadPeaks(const std::string& peakfile);

};

#endif  // SRC_PEAKINTERVALS_H__
