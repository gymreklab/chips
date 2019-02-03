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
  peakfiletype = "";
  countindex = -1;

  // Required learn inputs
  chipbam = "";
  remove_pct = 0;

  // Simulation experiment parameters
  numcopies = 100;
  numreads = 1000000;
  readlen = 36;
  paired = false;
  n_threads = 1;

  // Simulation model parameters
  gamma_alpha = 15.7;
  gamma_beta = 15.5;
  ratio_s = 0.001;
  ratio_f = 0.001;

  // Additional simulation parameters
  region = "";
  binsize = 100000;

  // sequencer
  sequencer_type = "";
  sub_rate = 0;
  ins_rate = 0;
  del_rate = 0;

  // pcr
  pcr_rate = 1.0;
}

Options::~Options() {}
