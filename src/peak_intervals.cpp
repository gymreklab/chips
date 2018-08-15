#include "src/common.h"
#include "src/peak_intervals.h"

PeakIntervals::PeakIntervals(const std::string peakfile, const std::string peakfileType,
        const std::string bamfile, const std::int32_t count_colidx) {
  if (!LoadPeaks(peakfile, peakfileType, bamfile, count_colidx)) {
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
bool PeakIntervals::LoadPeaks(const std::string peakfile, const std::string peakfileType,
        const std::string bamfile, const std::int32_t count_colidx) {
  PeakLoader peakloader(peakfile, peakfileType, bamfile, count_colidx);

  bool dataLoaded = peakloader.Load(peaks);

  if (dataLoaded){
    max_coverage = -1;
    for (int peakIndex=0; peakIndex<peaks.size(); peakIndex++){
      if (peaks[peakIndex].score > max_coverage){
        max_coverage = peaks[peakIndex].score;
      }
    }
  } 
  return dataLoaded;
}


/*
 * Inputs:
 *   - Fragment frag: a fragment
 *   - std::vector<Fragment>: a list of peaks
 *
 * Outputs:
 *   - float: probablity of the input fragment being bound
 * */
float PeakIntervals::SearchList(const Fragment& frag){
  std::vector<float> probBoundList;
  for(int peakIndex=0; peakIndex < peaks.size(); peakIndex++){
    if ((frag.chrom == peaks[peakIndex].chrom) &&\
            (frag.start >= peaks[peakIndex].start) &&\
            ( (frag.start+frag.length) <= (peaks[peakIndex].start+peaks[peakIndex].length))){
      probBoundList.push_back(peaks[peakIndex].score/max_coverage);
    }
  }

  // suppose probBoundList = [a, b, c]
  // result = 1- (1-a)*(1-b)*(1-c)
  if (probBoundList.size() > 0){
    float probBound = 1;
    for (int peakIndex=0; peakIndex < probBoundList.size(); peakIndex++){
       probBound *= (1-probBoundList[peakIndex]);
    }
    probBound = 1-probBound;
    return probBound;
  }else{
    return 0;
  }
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

  return SearchList(frag);

}
