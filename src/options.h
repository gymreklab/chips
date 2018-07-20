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

  // Required learn inputs
  std::string chipbam;

  // Simulation experiment parameters
  int numcopies;
  int numreads;
  int readlen;
  bool paired;

  // Simulation model parameters
  float gamma_alpha;
  float gamma_beta;

  virtual ~Options();
};

#endif  // SRC_OPTIONS_H__
