#include <cstring>
#include <iostream>
#include <stdlib.h>
#include <vector>

#include "src/bingenerator.h"
#include "src/common.h"
#include "src/fragment.h"
#include "src/library_constructor.h"
#include "src/options.h"
#include "src/pulldown.h"
#include "src/sequencer.h"
#include "src/stringops.h"

using namespace std;

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
       (PARAMETER_CHECK("--help", 6, parameterLength))) {
      showHelp = true;
    }
  }
  if (showHelp) {simulate_reads_help();}

  // do some parsing (all of these parameters require 2 strings)
  for (int i = 1; i < argc; i++) {
    int parameterLength = (int)strlen(argv[i]);

    if (PARAMETER_CHECK("-p", 2, parameterLength)) {
      if ((i+1) < argc) {
	options.peaksbed = argv[i+1];
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
    } else if (PARAMETER_CHECK("--numcopies", 11, parameterLength)) {
      if ((i+1) < argc) {
	options.numcopies = atoi(argv[i+1]);
	i++;
      }
    } else if (PARAMETER_CHECK("--numreads", 10, parameterLength)) {
      if ((i+1) < argc) {
	options.numreads = atoi(argv[i+1]);
	i++;
      }
    } else if (PARAMETER_CHECK("--readlen", 9, parameterLength)) {
      if ((i+1) < argc) {
	options.readlen = atoi(argv[i+1]);
	i++;
      }
    } else if (PARAMETER_CHECK("--gamma-frag", 12, parameterLength)) {
      if ((i+1) < argc) {
	std::vector<std::string> params;
	split_by_delim(argv[i+1], ',', params);
	options.gamma_alpha = atof(params[0].c_str());
	options.gamma_beta = atof(params[1].c_str());
	i++;
      }
    } else if (PARAMETER_CHECK("--spot", 6, parameterLength)) {
      if ((i+1) < argc) {
	options.ratio_s = atof(argv[i+1]);
	i++;
      }
    } else if (PARAMETER_CHECK("--frac", 6, parameterLength)) {
      if ((i+1) < argc) {
	options.ratio_f = atof(argv[i+1]);
	i++;
      }
    } else if (PARAMETER_CHECK("--region", 8, parameterLength)) {
      if ((i+1) < argc) {
	options.region = argv[i+1];
	i++;
      }
    } else if (PARAMETER_CHECK("--binsize", 9, parameterLength)) {
      if ((i+1) < argc) {
	options.binsize = atoi(argv[i+1]);
	i++;
      }
    } else if (PARAMETER_CHECK("--paired", 8, parameterLength)) {
      options.paired = true;
    } else {
      cerr << endl << "*****ERROR: Unrecognized parameter: " << argv[i] << " *****" << endl << endl;
      showHelp = true;
    }
  }

  // Check inputs
  if (options.peaksbed.empty()) {
    cerr << "****** ERROR: Must specify peaks with -p ******" << endl;
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
    vector<Fragment> pulldown_fragments, lib_fragments;

    // Perform in bins so we don't keep everything in memory at once
    BinGenerator bingenerator(options);

    while(bingenerator.GotoNextBin()) {
      /*** Step 1/2: Shearing + Pulldown ***/
      Pulldown pulldown(options, bingenerator.GetCurrentBin());
      pulldown.Perform(&pulldown_fragments);
    
      /*** Step 3: Library construction ***/
      LibraryConstructor lc(options);
      lc.Perform(pulldown_fragments, &lib_fragments);

      /*** Step 4: Sequencing ***/
      Sequencer seq(options);
      seq.Sequence(lib_fragments);
    }

    return 0;
    /******************************************************/
  } else {
    simulate_reads_help();
    return 1;
  }
}

void simulate_reads_help(void) {
  Options options;
  cerr << "\nTool:    asimon simreads" << endl;
  cerr << "Version: " << _GIT_VERSION << "\n";    
  cerr << "Summary: Simulate ChIP-seq reads for a set of peaks." << endl << endl;
  cerr << "Usage:   " << PROGRAM_NAME << " simreads -p peaks.bed -f ref.fa -o outprefix [OPTIONS] " << endl;
  cerr << "\n[Required arguments]: " << "\n";
  cerr << "     -p <peaks.bed>: BED file with peak regions" << "\n";
  cerr << "     -f <ref.fa>: FASTA file with reference genome" << "\n";
  cerr << "     -o <outprefix>: Prefix for output files" << "\n";
  cerr << "\n[Experiment parameters]: " << "\n";
  cerr << "     --numcopies <int>: Number of copies of the genome to simulate.\n"
       << "                        Default: " << options.numcopies << "\n";
  cerr << "     --numreads <int> : Number of reads (or read pairs) to simulate.\n"
       << "                        Default: " << options.numreads << "\n";
  cerr << "     --readlen <int>  : Read length to generate.\n"
       << "                        Default: " << options.readlen << "\n";
  cerr << "     --paired         : Simulate paired-end reads.\n"
       << "                        Default: false \n";
  cerr << "\n[Model parameters]: " << "\n";
  cerr << "     --gamma-frag <float>,<float>: Parameters for fragment length distribution (alpha, beta).\n"
       << "                                   Default: " << options.gamma_alpha << ","
       << options.gamma_beta << "\n";
  cerr << "     --spot <float>              : SPOT score (fraction of reads in peaks) \n"
       << "                                   Default: " << options.ratio_s << "\n";
  cerr << "     --frac <float>              : Fraction of the genome that is bound \n"
       << "                                   Default: " << options.ratio_f << "\n";
  cerr << "\n[Other options]: " << "\n";
  cerr << "     --region <str>              : Only simulate reads from this region chrom:start-end\n"
       << "                                   Default: genome-wide \n";
  cerr << "     --binsize <int>             : Consider bins of this size when simulating\n"
       << "                                 : Default: " << options.binsize << "\n";
  cerr << "\n";
  exit(1);
}
