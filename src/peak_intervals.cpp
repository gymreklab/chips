#include "src/common.h"
#include "src/peak_intervals.h"

PeakIntervals::PeakIntervals(const std::string& peakfile) {
  if (!LoadPeaks(peakfile)) {
    PrintMessageDieOnError("Error loading peaks from " + peakfile, M_ERROR);
  }
}

PeakIntervals::~PeakIntervals() {}

/*
  Inputs:
  - string peakfile: filename of peaks. tab-separated with chrom, start, end, score

  Outputs:
  - bool: return true if successful, false if error loading peaks

  This function puts the peaks into a searachable data structure
 */
bool PeakIntervals::LoadPeaks(const std::string& peakfile) {
  // TODO An/Michael
  return false;
}

/*
  Inputs:
  - Fragment frag: has chrom, start, and length of a fragment

  Outputs:
  - float: probability that the fragment is bound

  If the fragment doesn't overlap a peak, return 0
  If it overlaps one ore more peak, return the max score across all peaks
 */
float PeakIntervals::GetOverlap(const Fragment& frag) {
  // TODO An/Michael
  return 0;
}
