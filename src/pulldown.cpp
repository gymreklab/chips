#include "src/pulldown.h"

#include <iostream>
#include <random>

int Pulldown::peakIndexStart = 0;
std::string Pulldown::prev_chrom = "";

Pulldown::Pulldown(const Options& options, const GenomeBin& gbin) {
  chrom = gbin.chrom;
  start = gbin.start;
  end = gbin.end;
  numcopies = options.numcopies;
  gamma_alpha = options.gamma_alpha;
  gamma_beta = options.gamma_beta;
  ratio_beta = options.ratio_f*(1-options.ratio_s)/(options.ratio_s*(1-options.ratio_f));

  debug_pulldown = true; // TODO remove
  if (debug_pulldown) {
    PrintMessageDieOnError("Loading peaks", M_DEBUG);
  }
}

void Pulldown::Perform(vector<Fragment>* output_fragments, PeakIntervals* pintervals) {
  // Set up 
  std::default_random_engine generator;
  std::gamma_distribution<float> fragdist(gamma_alpha, gamma_beta);
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
  for (int i = 0; i < numcopies; i++) {
    pintervals->resetSearchScope(peakIndexStart); 
    if (debug_pulldown) {
      PrintMessageDieOnError("Genome copy " + to_string(i), M_DEBUG);
    }
    current_pos = start;
    // Break up into fragment lengths drawn from gamma distribution
    while (current_pos < end) {
      fsize = fragdist(generator);
      fstart = current_pos; fend = current_pos+fsize;
      if (fend > end) {break;}
      Fragment frag(chrom, current_pos, fsize);
      peak_score = pintervals->GetOverlap(frag);
      bound = (rand()/double(RAND_MAX) < peak_score);
      // TODO debug below
      if (peak_score > 0 && debug_pulldown) {
	cerr << chrom << " " << fstart << " " << fend << " " << " " << peak_score << " " << bound << " " << ratio_beta << endl;
      }
      if (bound) {
	output_fragments->push_back(frag); // alpha=1
      } else {
	if (rand()/double(RAND_MAX) < (ratio_beta * (pintervals->prob_pd_given_b) )) {
	  output_fragments->push_back(frag);
	}
      }
      current_pos += fsize;
    }
  }

  peakIndexStart = pintervals->peakIndexStart;
}

