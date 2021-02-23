#include "options.h"

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
  noscale = false;
  scale_outliers = false;

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
  gamma_k = 15.67;
  gamma_theta = 15.49;
  ratio_s = 0.17594;
  ratio_f = 0.03713;
  pcr_rate = 1.0;
  model_file = "";
  recompute_f = false;

  // random seed
  seed = 0;

  // Additional simulation parameters
  region = "";
  binsize = 100000;

  // Additional learn parameters
  downsample = 1.0;
  skip_frag = false;
  skip_pd = false;
  skip_pcr = false;
  output_frag_lens = false;

  // sequencer
  sequencer_type = "";
  sub_rate = 0;
  ins_rate = 0;
  del_rate = 0;

  // learn only: intensity threshold
  intensity_threshold = 100;
  intensity_threshold_scale = 0.0;
  estimate_frag_length = 300;

  // Other options
  verbose = false;
}

Options::~Options() {}
