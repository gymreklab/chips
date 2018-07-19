#include <cstring>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <sstream>

#include "src/bam_io.h"
#include "src/common.h"
#include "src/options.h"

using namespace std;

// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)
const bool DEBUG = true;

// Function declarations
void learn_help(void);
bool learn_frag(const std::string& bamfile, float* param);

bool learn_frag(const std::string& bamfile, float* param) {
  // TODO which params do we need to learn? probably more than one  
  /*
    Learn fragment length distribution from an input BAM file
    Fragment lengths follow a gamma distribution.

    Inputs:
    - bamfile (std::string): path to the BAM file
    Outputs:
    - param (float): parameter of gamma distribution (TODO update with correct name and number of params to be learned)
   */

  /* First, get a vector of the fragment lengths */
  int maxreads = 10000; int numreads = 0; // don't look at more than this many reads
  BamCramReader bamreader(bamfile);
  const BamHeader* bamheader = bamreader.bam_header();
  // Get first chrom to look at fragment lengths
  std::vector<std::string> seq_names = bamheader->seq_names();
  std::vector<uint32_t> seq_lengths = bamheader->seq_lengths();
  if (seq_names.size() > 0 && seq_lengths.size() > 0) {
    bamreader.SetRegion(seq_names[0], 0, seq_lengths[0]);
  } else {
    return false;
  }
  std::vector<int32_t> fraglengths;
  int32_t tlen;
  BamAlignment aln;
  while(bamreader.GetNextAlignment(aln) && numreads<maxreads) {
    tlen = aln.TemplateLength();
    if (tlen > 0) {
      fraglengths.push_back(abs(tlen));
      numreads++;
    }
    //    cerr << abs(tlen) << endl; // if you want to print out for debugging
  }

  /* Now, fit fraglengths to a gamma distribution */
  // TODO (MICHAEL/AN) FILL THIS PART IN, RETURN TRUE IF SUCCESSFUL
  // MAY NEED TO CHANGE fraglenths DATA STRUCTURE depending on what optimizer takes
  *param = 0; // example for how to set result
  if (DEBUG) {
    std::stringstream ss;
    ss << "Learned fragment length param " << *param;
    PrintMessageDieOnError(ss.str(), M_DEBUG);
  }
  return true;
}

int learn_main(int argc, char* argv[]) {
  bool showHelp = false;
  Options options;

  // check to see if we should print out some help
  if(argc <= 1) showHelp = true;
  for(int i = 1; i < argc; i++) {
    int parameterLength = (int)strlen(argv[i]);

    if ((PARAMETER_CHECK("-h", 2, parameterLength)) ||
       (PARAMETER_CHECK("--help", 5, parameterLength))) {
      showHelp = true;
    }
  }
  if (showHelp) {learn_help();}

  // do some parsing (all of these parameters require 2 strings)
  for (int i = 1; i<argc; i++) {
    int parameterLength = (int)strlen(argv[i]);
    if (PARAMETER_CHECK("-b", 2, parameterLength)) {
      if ((i+1) < argc) {
	options.chipbam = argv[i+1];
	i++;
      }
    } else if (PARAMETER_CHECK("-p", 2, parameterLength)) {
      if ((i+1) < argc) {
	options.peaksbed = argv[i+1];
	i++;
      }
    } else if (PARAMETER_CHECK("-o", 2, parameterLength)) {
      if ((i+1) < argc) {
	options.outprefix = argv[i+1];
	i++;
      }
    } else {
      cerr << endl << "******ERROR: Unrecognized parameter: " << argv[i] << " ******" << endl << endl;
      showHelp = true;
    }
  }

  // Check inputs
  if (options.chipbam.empty()) {
    cerr << "****** ERROR: Must specify BAM file with -b ******" << endl;
    showHelp = true;
  }
  if (options.peaksbed.empty()) {
    cerr << "****** ERROR: Must specify peaks file with -p ******" << endl;
    showHelp = true;
  }
  if (options.outprefix.empty()) {
    cerr << "****** ERROR: Must specify outprefix with -o ******" << endl;
    showHelp = true;
  }

  if (!showHelp) {
    /***************** Main implementation ***************/

    /*** Learn fragment size disbribution parameters ***/
    float frag_param;
    if (!learn_frag(options.chipbam, &frag_param)) {
      PrintMessageDieOnError("Error learning fragment length distribution", M_ERROR);
    }

    /*** Learn pulldown ratio parameters ***/
    // TODO

    /*** Write params to file ***/
    // TODO

    return 0;
  } else {
    learn_help();
    return 1;
  }
}

void learn_help(void) {
  cerr << "\nTool:    asimon learn" << endl;
  cerr << "Version: " << _GIT_VERSION << "\n";    
  cerr << "Summary: Learn parameters from a ChIP dataset." << endl << endl;
  cerr << "Usage:   " << PROGRAM_NAME << " learn -b reads.bam -p peak.bed -o outprefix [OPTIONS] " << endl << endl;
  cerr << "[Required arguments]: " << "\n";
  cerr << "         -b <reads.bam>: BAM file with ChIP reads (.bai index required)" << "\n";
  cerr << "         -p <peaks.bed>: BED file with peak regions (Homer format)" << "\n";
  cerr << "         -o <outprefix>: Prefix for output files" << "\n";
  cerr << "\n";
  exit(1);
}
