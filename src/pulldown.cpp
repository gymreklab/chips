#include "src/pulldown.h"
#include <chrono>

#include <iostream>
#include <random>

Pulldown::Pulldown(const Options& options, const GenomeBin& gbin,\
        std::string& _prev_chrom, int& _peakIndexStart, int& _start_offset) {
  chrom = gbin.chrom;
  start = gbin.start;
  end = gbin.end;
  numcopies = options.numcopies;
  gamma_k = options.gamma_k;
  gamma_theta = options.gamma_theta;
  ratio_beta = options.ratio_f*(1-options.ratio_s)/(options.ratio_s*(1-options.ratio_f));

  prev_chrom = _prev_chrom;
  peakIndexStart = _peakIndexStart;
  start_offset_ptr = & _start_offset;
}

void Pulldown::Perform(vector<Fragment>* output_fragments, PeakIntervals* pintervals) {
  // Set up
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator(seed);
  srand(seed);
  std::gamma_distribution<float> fragdist(gamma_k, gamma_theta);
  int32_t current_pos;
  int32_t fstart, fend;
  float fsize;
  bool bound;
  float peak_score;

  // update the start index of peaks
  if (chrom != prev_chrom){
    peakIndexStart = 0;
    prev_chrom = chrom;
  }

  // Perform separate shearing for each copy of the genome
  // pintervals->resetSearchScope(peakIndexStart); 
  int peakIndex = peakIndexStart;
  current_pos = start + *start_offset_ptr;
  // Break up into fragment lengths drawn from gamma distribution
  while (current_pos < end) {
    fsize = fragdist(generator);
    fstart = current_pos; fend = current_pos+fsize;
    if (fend > end) {
      if (fstart < end){
        *start_offset_ptr = fend-end;
      }else{
        break;
      }
    }
    Fragment frag(chrom, current_pos, fsize);
    peak_score = pintervals->GetOverlap(frag, peakIndex);

    bound = (rand()/double(RAND_MAX) < peak_score);
    if (bound) {
      output_fragments->push_back(frag);
    } else{
      if (rand()/double(RAND_MAX) < ratio_beta) {
          output_fragments->push_back(frag);
      }
    }
    current_pos += fsize;
  }

  peakIndexStart = peakIndex;
}

