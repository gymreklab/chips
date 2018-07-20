#include "src/options.h"

using namespace std;

/*
  Initialize options
 */
Options::Options() {
  // Required simulation inputs
  reffa = "";
  peaksbed = "";
  outprefix = "";

  // Required learn inputs
  chipbam = "";

  // Simulation experiment parameters
  numcopies = 100;
  numreads = 1000000;
  readlen = 36;
  paired = false;

  // Simulation model parameters
  gamma_alpha = 100; // TODO set to more reasonable defaults
  gamma_beta = 0.5;
}

Options::~Options() {}
