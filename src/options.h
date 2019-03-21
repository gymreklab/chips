#ifndef SRC_OPTIONS_H__
#define SRC_OPTIONS_H__

#include <vector>
#include <string>

#include <stdint.h>

class Options {
 public:
  Options();

  // Required simulation inputs
  std::string reffa;
  std::string peaksbed;
  std::string outprefix;
  std::string peakfiletype;
  std::int32_t countindex;
  bool noscale;
  bool scale_outliers;

  // Required learn inputs
  std::string chipbam;
  float remove_pct;

  // Simulation experiment parameters
  int numcopies;
  int numreads;
  int readlen;
  bool paired;
  int n_threads;

  // Simulation model parameters
  float gamma_k;
  float gamma_theta;
  float ratio_s;
  float ratio_f;
  float pcr_rate; 
  std::string model_file;
  bool recompute_f;

  std::string sequencer_type;
  float sub_rate;
  float del_rate;
  float ins_rate;
  // Additional simulation parameters
  std::string region;
  int binsize;

  int intensity_threshold;
  int estimate_frag_length;

  // Learn options
  bool skip_frag;
  bool output_frag_lens;
  virtual ~Options();

  // Other options
  bool verbose;
};

#endif  // SRC_OPTIONS_H__
