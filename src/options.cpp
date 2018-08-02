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
  gamma_alpha = 100;
  gamma_beta = 0.5;
  ratio_s = 0.001;
  ratio_f = 0.00001;

  // Additional simulation parameters
  region = "";
  binsize = 100000;
}

Options::~Options() {}
