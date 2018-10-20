#include <cstring>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <sstream>
#include <math.h>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/polygamma.hpp>

#include "src/bam_io.h"
#include "src/common.h"
#include "src/options.h"
#include "src/peak_loader.h"
using namespace std;

// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)
const bool DEBUG = true;

// Function declarations
void learn_help(void);
bool learn_ratio(const std::string& bamfile, const std::string& peakfile,
                    const std::string& peakfileType, const std::int32_t count_colidx,
                    const float remove_pct, float* ab_ratio_ptr, 
                    float *s_ptr, float* f_ptr);
bool compare_location(Fragment a, Fragment b);
bool learn_frag(const std::string& bamfile, float* alpha, float* beta);

bool learn_frag(const std::string& bamfile, float* alpha, float* beta) {
  /*
    Learn fragment length distribution from an input BAM file
    Fragment lengths follow a gamma distribution.

    Inputs:
    - bamfile (std::string): path to the BAM file
    Outputs:
    - alpha (float): parameter of gamma distribution
    - beta  (float): parameter of gamma distribution
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

  /* Now, fit fraglengths to a gamma distribution.
     Use Maximum Likelihood Estimation to estimate the
     Gamma Distribution parameters */
 
  const float EPSILON = 1e-4; // Value to check for convergence
  
  float total_frag_len = 0;     // sum of all the frag lengths
  float total_log = 0;          // sum of log of each frag length

  // get sum of all frag lengths and log sum of each frag length
  for (int frag = 0; frag < fraglengths.size(); frag++)
  {
    total_frag_len += fraglengths[frag];
    total_log += log(fraglengths[frag]);
  }

  // mean of frag lengths and log mean of fraglengths
  float mean_frag_length = total_frag_len/fraglengths.size();
  float total_log_mean = total_log/fraglengths.size();
  
  // Starting point for the value of a
  float a = 0.5/(log(mean_frag_length) - total_log_mean);
  float new_a = 0;
  
  // estimate the value for a using maximum likelihood estimate
  while (true)
  {
    // evaluate updated a
    float update = (1/a) + ((total_log_mean - log(mean_frag_length) + log(a)
                 - boost::math::digamma(a))/(a - a*a*boost::math::polygamma(1, a)));
    new_a = 1 / update;

    // a converges
    if (abs(new_a - a) < EPSILON)
      break;

    a = new_a;
  }

  *beta = mean_frag_length/a;
  *alpha = a;
  
  if (DEBUG) {
    std::stringstream ss;
    ss << "Learned fragment length params alpha: " << *alpha << " and beta: " << *beta;
    PrintMessageDieOnError(ss.str(), M_DEBUG);
  }
  return true;
}


bool learn_ratio(const std::string& bamfile, const std::string& peakfile,
        const std::string& peakfileType, const std::int32_t count_colidx, 
        const float remove_pct, float* ab_ratio_ptr, float *s_ptr, float* f_ptr){
  /*
    Learn the ratio of alpha to beta from an input BAM file,
    and its correspounding peak file.

    Here,
    alpha = Prob(P|B); beta = Prob(P| NOT B)
    f:  fraction of the genome bound,
        approximated as the length of all peaks divided by length of the genome
    s:  SPOT score, or fraction of all reads falling peaks
    
    ab_ratio = (s/(1-s)) * ((1-f)/f)

    Inputs:
    - bamfile (std::string): path to the BAM file
    - peakfile (std::string): path to the ChIP-Seq peak file
    Outputs:
    - ab_ratio (float): the ratio of alpha to beta
   */


  // Read peak locations from the ChIP-seq file,
  // and calculate the total length of peaks across the genomes
  std::vector<Fragment> peaks;
  PeakLoader peakloader(peakfile, peakfileType, bamfile, count_colidx);
  peakloader.Load(peaks);

  // Remove top remove_pct% of peaks default is do not remove
  if (remove_pct > 0)
  {
    int keep_peaks = floor(peaks.size()*(1 - remove_pct));
    auto score_min = [](Fragment a, Fragment b) {return a.score < b.score;};
    std::sort(peaks.begin(), peaks.end(), score_min);
    peaks.resize(keep_peaks, Fragment("", 0, 0));
  }

  int plen = 0;
  cout<<peaks.size()<<endl;
  for(int peak_index = 0; peak_index < peaks.size(); peak_index++){
    plen += peaks[peak_index].length;
  }

  // calculate f and s, and then ab_ratio
  float f = (float)plen / (float) (peakloader.total_genome_length);
  float s = (float)(peakloader.tagcount_in_peaks) / (float)(peakloader.total_tagcount);
  std::cout <<"peak-length: "<<plen<<"\ttotal: "<< peakloader.total_genome_length<<std::endl;
  std::cout <<"#reads in peaks: "<<peakloader.tagcount_in_peaks<<"\t#reads: "<<peakloader.total_tagcount<<std::endl;

  // outputs
  float ab_ratio = (s/(1-s)) * ((1-f)/f);
  *ab_ratio_ptr = ab_ratio;
  *f_ptr = f;
  *s_ptr = s;
  return true;
}

bool compare_location(Fragment a, Fragment b){
  if (a.chrom != b.chrom){
    return (a.chrom < b.chrom);
  }else if (a.start != b.start){
    return (a.start < b.start);
  }else if (a.length != b.length){
    return (a.length < b.length);
  }else{
    return true;
  }
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
    } else if (PARAMETER_CHECK("-t", 2, parameterLength)){
      if ((i+1) < argc) {
    options.peakfiletype = argv[i+1];
	i++;
      }
    } else if (PARAMETER_CHECK("-c", 2, parameterLength)){
      if ((i+1) < argc) {
    options.countindex = std::atoi(argv[i+1]);
	i++;
      }
    } else if (PARAMETER_CHECK("-r", 2, parameterLength)){
      if ((i+1) < argc) {
    options.remove_pct = (float) (std::atof(argv[i+1])/100);
    i++;
      }
    } else {
      cerr << endl << "******ERROR: Unrecognized parameter: " << argv[i] << " ******" << endl << endl;
      showHelp = true;
    }
  }

  // Check inputs
  if (options.peaksbed.empty()) {
    cerr << "****** ERROR: Must specify peaks file with -p ******" << endl;
    showHelp = true;
  }
  if (options.outprefix.empty()) {
    cerr << "****** ERROR: Must specify outprefix with -o ******" << endl;
    showHelp = true;
  }
  if (options.peakfiletype.empty()) {
    cerr << "****** ERROR: Must specify peakfiletype with -t ******" << endl;
    showHelp = true;
  }


  if (!showHelp) {
    /***************** Main implementation ***************/

    /*** Learn fragment size disbribution parameters ***/
    float frag_param_a;
    float frag_param_b;
    if (!learn_frag(options.chipbam, &frag_param_a, &frag_param_b)) {
      PrintMessageDieOnError("Error learning fragment length distribution", M_ERROR);
    }

    /*** Learn pulldown ratio parameters ***/
    float ab_ratio;
    float s;
    float f;
    if (!learn_ratio(options.chipbam, options.peaksbed, options.peakfiletype, options.countindex, options.remove_pct,
        &ab_ratio, &s, &f)){
      PrintMessageDieOnError("Error learning pulldown ratio", M_ERROR);
    }
    cout << "ab_ratio: " << ab_ratio << endl;
    cout << "f: " << f << endl;
    cout << "s: " << s << endl;

    /*** Write params to file ***/
    // TODO

    return 0;
  } else {
    learn_help();
    return 1;
  }
}

//TODO add in -r argument not required
void learn_help(void) {
  cerr << "\nTool:    asimon learn" << endl;
  cerr << "Version: " << _GIT_VERSION << "\n";    
  cerr << "Summary: Learn parameters from a ChIP dataset." << endl << endl;
  cerr << "Usage:   " << PROGRAM_NAME << " learn -b reads.bam -p peak.bed -o outprefix [OPTIONS] " << endl << endl;
  cerr << "[Required arguments]: " << "\n";
  cerr << "         -b <reads.bam>: BAM file with ChIP reads (.bai index required)" << "\n";
  cerr << "         -p <peaks.bed>: BED file with peak regions (Homer format)" << "\n";
  cerr << "         -t <peakfile_type>: File type of the input peak file" << "\n";
  cerr << "         -o <outprefix>: Prefix for output files" << "\n";
  cerr << "\n";
  exit(1);
}
