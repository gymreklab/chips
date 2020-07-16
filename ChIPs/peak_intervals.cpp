#include "common.h"
#include "peak_intervals.h"
#include "ref_genome.h"

PeakIntervals::PeakIntervals(const Options& options, const std::string peakfile, const std::string peakfileType,
			     const std::string bamfile, const std::int32_t count_colidx) {
  if (!LoadPeaks(options, peakfile, peakfileType, bamfile, count_colidx)) {
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
bool PeakIntervals::LoadPeaks(const Options& options,
			      const std::string peakfile, const std::string peakfileType,
			      const std::string bamfile, const std::int32_t count_colidx) {
  PeakLoader peakloader(peakfile, peakfileType, bamfile, count_colidx);

  std::vector<Fragment> peaks;
  float frag_length = options.gamma_k * options.gamma_theta; 
  bool dataLoaded = peakloader.Load(peaks, options.region, frag_length, options.noscale, options.scale_outliers);

  if (dataLoaded){
    total_bound_length = 0;
    for (int peakIndex=0; peakIndex<peaks.size(); peakIndex++){
      total_bound_length += (peaks[peakIndex].length * peaks[peakIndex].score);
    }
    total_genome_length = peakloader.total_genome_length;

    // convert the vector into a map with keys as chromosome names
    for (int peakIndex=0; peakIndex<peaks.size(); peakIndex++){
      peak_map[peaks[peakIndex].chrom].push_back(peaks[peakIndex]);
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
float PeakIntervals::SearchList(const Fragment& frag, int& peakIndexStart){
  std::vector<float> probBoundList;
  std::int32_t frag_start, frag_end, peak_start, peak_end;
  float overlap;
  std::vector<Fragment> & peaks = peak_map[frag.chrom];
  for(int peakIndex=peakIndexStart; peakIndex < peaks.size(); peakIndex++){
    if (frag.chrom == peaks[peakIndex].chrom){
      frag_start = frag.start;
      frag_end = frag.start+frag.length;
      peak_start = peaks[peakIndex].start;
      peak_end = peaks[peakIndex].start+peaks[peakIndex].length;

      if (frag_end <= peak_start){
        break; // there remain no overlapped peaks
      }else if(frag_start >= peak_end){
        peakIndexStart += 1; // move to the next peak
      }else{
        overlap = (float) (std::min(peak_end,frag_end) - std::max(peak_start, frag_start)) / (float)(frag_end-frag_start);
        probBoundList.push_back(overlap*peaks[peakIndex].score);
      }
    }else{
        std::cerr << "****** ERROR: Unexpected errors in PeakInterval/SearchList ******" << std::endl;
        std::exit(1);
    }
  }
  
  // suppose probBoundList = [a, b, c]
  // result = 1- (1-a)*(1-b)*(1-c)
  if (probBoundList.size() > 0){
    //std::cout<<"found"<<std::endl;
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
float PeakIntervals::GetOverlap(const Fragment& frag, int& peakIndexStart) {
  float score = SearchList(frag, peakIndexStart);
  return score;
}

