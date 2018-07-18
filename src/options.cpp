#include "src/options.h"

using namespace std;

/*
  Initialize options
 */
Options::Options() {
  // Required simulation inputs
  reffa = "";
  regionsbed = "";
  outprefix = "";

  // Required learn inputs
  chipbam = "";
}

Options::~Options() {}
