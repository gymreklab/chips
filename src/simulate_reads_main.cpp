#include <cstring>
#include <iostream>
#include <stdlib.h>
#include <vector>

#include "src/fragment.h"
#include "src/library_constructor.h"
#include "src/options.h"
#include "src/pulldown.h"
#include "src/sequencer.h"
#include "src/shearer.h"

using namespace std;

// define our program name
#define PROGRAM_NAME "asimon simreads"

// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

// Function declarations
void simulate_reads_help(void);

int simulate_reads_main(int argc, char* argv[]) {
  bool showHelp = false;
  Options options;

  // check to see if we should print out some help
  if(argc <= 1) showHelp = true;
  for (int i = 1; i < argc; i++) {
    int parameterLength = (int)strlen(argv[i]);

    if ((PARAMETER_CHECK("-h", 2, parameterLength)) ||
       (PARAMETER_CHECK("--help", 5, parameterLength))) {
      showHelp = true;
    }
  }
  if (showHelp) {simulate_reads_help();}

  // do some parsing (all of these parameters require 2 strings)
  for (int i = 1; i < argc; i++) {
    int parameterLength = (int)strlen(argv[i]);

    if (PARAMETER_CHECK("-r", 2, parameterLength)) {
      if ((i+1) < argc) {
	options.regionsbed = argv[i+1];
	i++;
      }
    } else if (PARAMETER_CHECK("-f", 2, parameterLength)) {
      if ((i+1) < argc) {
	options.reffa = argv[i+1];
	i++;
      }
    } else if (PARAMETER_CHECK("-o", 2, parameterLength)) {
      if ((i+1) < argc) {
	options.outprefix = argv[i+1];
	i++;
      }
    } else {
      cerr << endl << "*****ERROR: Unrecognized parameter: " << argv[i] << " *****" << endl << endl;
      showHelp = true;
    }
  }

  // Check inputs
  if (options.regionsbed.empty()) {
    cerr << "****** ERROR: Must specify regions with -r ******" << endl;
    showHelp = true;
  }
  if (options.reffa.empty()) {
    cerr << "****** ERROR: Must specify reffa with -f ******" << endl;
    showHelp = true;    
  }
  if (options.outprefix.empty()) {
    cerr << "****** ERROR: Must specify outprefix with -o ******" << endl;
    showHelp = true;
  }

  if (!showHelp) {
    /***************** Main implementation ***************/
    // Set up
    vector<Fragment> sheared_fragments, pulldown_fragments, lib_fragments;

    /*** Step 1: Shearing ***/
    Shearer shearer(options);
    shearer.Shear(&sheared_fragments);

    /*** Step 2: Pulldown ***/
    Pulldown pulldown(options);
    pulldown.Perform(sheared_fragments, &pulldown_fragments);
    
    /*** Step 3: Library construction ***/
    LibraryConstructor lc(options);
    lc.Perform(pulldown_fragments, &lib_fragments);

    /*** Step 4: Sequencing ***/
    Sequencer seq(options);
    seq.Sequence(lib_fragments);

    return 0;
    /******************************************************/
  } else {
    simulate_reads_help();
    return 1;
  }
}

void simulate_reads_help(void) {
  cerr << "\nTool:    asimon simreads" << endl;
  cerr << "Version: " << _GIT_VERSION << "\n";    
  cerr << "Summary: Simulate ChIP-seq reads for a set of peaks." << endl << endl;
  cerr << "Usage:   " << PROGRAM_NAME << " -r regions.bed -f ref.fa -o outprefix [OPTIONS] " << endl << endl;
  cerr << "[Required arguments]: " << "\n";
  cerr << "     -r <regions.bed>: BED file with peak regions" << "\n";
  cerr << "     -f <ref.fa>: FASTA file with reference genome" << "\n";
  cerr << "     -o <outprefix>: Prefix for output files" << "\n";
  cerr << "\n";
  exit(1);
}
