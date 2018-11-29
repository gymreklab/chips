#include "src/common.h"
#include "src/peak_intervals.h"
#include "src/ref_genome.h"

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
  bool dataLoaded = peakloader.Load(peaks, options.region);

  if (dataLoaded){
    // calculate the maximum coverage
    max_coverage = 0;
    for (int peakIndex=0; peakIndex<peaks.size(); peakIndex++){
      if (peaks[peakIndex].score > max_coverage){
        max_coverage = peaks[peakIndex].score;
      }
    }

    // calculate Prob(pulled down|bound)
    double total_signals = 0;
    double total_signals_max = 0;
    for (int peakIndex=0; peakIndex<peaks.size(); peakIndex++){
      total_signals_max += (peaks[peakIndex].length * max_coverage);
      total_signals += (peaks[peakIndex].length * peaks[peakIndex].score);
    }
    prob_pd_given_b = total_signals/total_signals_max;

    // convert the vector into a map with keys as chromosome names
    for (int peakIndex=0; peakIndex<peaks.size(); peakIndex++){
      peak_map[peaks[peakIndex].chrom].push_back(peaks[peakIndex]);
    }
  }

  EstNumFrags(options, peaks);
  return dataLoaded;
}

void PeakIntervals::EstNumFrags(const Options& options, std::vector<Fragment> peaks){
    int total_length = 0;
    int length_b = 0;
    double numfrags_b = 0;
    double numfrags_ub = 0;
    double ratio_beta = options.ratio_f*(1-options.ratio_s)/(options.ratio_s*(1-options.ratio_f));
    double prob_pd_given_ub = ratio_beta * prob_pd_given_b;
    double frag_length = options.gamma_alpha * options.gamma_beta;
    //if(options.region.empty()){

    // estimate number of reads per run
    // bound regions
    for (int peakIndex=0; peakIndex<peaks.size(); peakIndex++){
      length_b += peaks[peakIndex].length;
      numfrags_b += ((double) peaks[peakIndex].length / frag_length) * (peaks[peakIndex].score/max_coverage);
    }

    // unbound regions
    if(options.region.empty()){
      std::map<std::string, int> chromLengths;
      RefGenome ref (options.reffa);
      if (!ref.GetLengths(&chromLengths))
        PrintMessageDieOnError("Could not gather chromosome lengths from "
                                + options.reffa, M_ERROR);
      for (auto _chrom : chromLengths){total_length += _chrom.second;}
    }else{
      // parse "chrID:start-end"
      std::vector<std::string> parts;
      std::stringstream ss(options.region);
      std::string split;
      while(std::getline(ss, split, ':')) parts.push_back(split);
      if (parts.size() != 2) PrintMessageDieOnError("Improper region input format should be chrom:start-end", M_ERROR);
      std::string region_chrom = parts[0];
      std::string start_end = parts[1];

      parts.clear();
      ss.str(start_end);
      ss.clear();

      while(getline(ss, split, '-')) parts.push_back(split);
      if (parts.size() != 2) PrintMessageDieOnError("Improper region input format should be chrom:start-end", M_ERROR);
      int region_start = stoi(parts[0]);
      int region_end = stoi(parts[1]);
      total_length = region_end - region_start;
    }
    numfrags_ub = prob_pd_given_ub * (double) (total_length - length_b) / frag_length;
    // put together
    double numfrags_per_run = numfrags_ub + numfrags_b;
    prob_frag_kept = (double) (options.numreads) / (double) (numfrags_per_run * options.numcopies);
    //std::cout<<numfrags_per_run<<"****"<<prob_frag_kept<<std::endl;
}

void PeakIntervals::resetSearchScope(const int index){
  peakIndexStart = index;
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
  uint32_t frag_start, frag_end, peak_start, peak_end;
  float overlap;
  /*
  for(int peakIndex=0; peakIndex < peaks.size(); peakIndex++){
    if (frag.chrom == peaks[peakIndex].chrom){
      frag_start = frag.start;
      frag_end = frag.start+frag.length;
      peak_start = peaks[peakIndex].start;
      peak_end = peaks[peakIndex].start+peaks[peakIndex].length;
      if (peak_start < frag_start){
        if (peak_end > frag_start){
            if (peak_end > frag_end){
              probBoundList.push_back(peaks[peakIndex].score/max_coverage);
            }else{
              overlap = (float)(peak_end-frag_start) / (float)(peak_end-peak_start);
              probBoundList.push_back(overlap*peaks[peakIndex].score/max_coverage);
            }
        }else{
            overlap = 0;
        }
      }else if (peak_start < frag_end){
        overlap = (float)( std::min(peak_end,frag_end) - peak_start)/(float)(peak_end-peak_start);
        probBoundList.push_back(overlap*peaks[peakIndex].score/max_coverage);
      }else{
        overlap = 0;
      }
    }
  }*/
  std::vector<Fragment> & peaks = peak_map[frag.chrom];
  for(int peakIndex=peakIndexStart; peakIndex < peaks.size(); peakIndex++){
    if (frag.chrom == peaks[peakIndex].chrom){
      frag_start = frag.start;
      frag_end = frag.start+frag.length;
      peak_start = peaks[peakIndex].start;
      peak_end = peaks[peakIndex].start+peaks[peakIndex].length;
      if (peak_start < frag_start){
        if (peak_end > frag_start){
            if (peak_end > frag_end){
              probBoundList.push_back(peaks[peakIndex].score/max_coverage);
            }else{
              overlap = (float)(peak_end-frag_start) / (float)(frag_end-frag_start);
              //overlap = 1;
              probBoundList.push_back(overlap*peaks[peakIndex].score/max_coverage);
            }
        }else{
            overlap = 0;
            peakIndexStart += 1;
        }
      }else if (peak_start < frag_end){
        overlap = (float)( std::min(peak_end,frag_end) - peak_start)/(float)(frag_end-frag_start);
        //overlap = 1;
        probBoundList.push_back(overlap*peaks[peakIndex].score/max_coverage);
      }else{
        overlap = 0;
        break;
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
float PeakIntervals::GetOverlap(const Fragment& frag) {
  float score = SearchList(frag);
  return score;
  
}
