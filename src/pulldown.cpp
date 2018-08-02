#include "src/pulldown.h"

#include <iostream>
#include <random>

Pulldown::Pulldown(const Options& options, const GenomeBin& gbin) {
  chrom = gbin.chrom;
  start = gbin.start;
  end = gbin.end;
  numcopies = options.numcopies;
  gamma_alpha = options.gamma_alpha;
  gamma_beta = options.gamma_beta;
  ratio_beta = options.ratio_f*(1-options.ratio_s)/(options.ratio_s*(1-options.ratio_f));
}

void Pulldown::Perform(vector<Fragment>* output_fragments) {
  // Set up 
  std::default_random_engine generator;
  std::gamma_distribution<float> fragdist(gamma_alpha, gamma_beta);
  int32_t current_pos;
  int32_t fstart, fend;
  float fsize;
  bool bound;

  // Perform separate shearing for each copy of the genome
  for (int i = 0; i < numcopies; i++) {
    current_pos = start;
    // Break up into fragment lengths drawn from gamma distribution
    while (current_pos < end) {
      fsize = fragdist(generator);
      fstart = current_pos; fend = current_pos+fsize;
      if (fend > end) {break;}
      Fragment frag(chrom, current_pos, fsize);
      bound = false; // TODO check for overlap peak and use score
      if (bound) {
	output_fragments->push_back(frag); // alpha=1
      } else {
	if (rand()/double(RAND_MAX) < ratio_beta) {
	  output_fragments->push_back(frag);
	}
      }
      current_pos += fsize;
    }
  }
}


Pulldown::~Pulldown() {}
