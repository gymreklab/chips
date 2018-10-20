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
  float gamma_alpha;
  float gamma_beta;
  float ratio_s;
  float ratio_f;

  // Additional simulation parameters
  std::string region;
  int binsize;

  virtual ~Options();
};

#endif  // SRC_OPTIONS_H__
