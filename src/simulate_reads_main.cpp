#include <cstring>
#include <iostream>
#include <fstream>
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
#include "src/peak_intervals.h"
#include "src/multithread.h"
#include "src/multithread.cpp"
const bool DEBUG_SIM=true;

// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

// Function declarations
void simulate_reads_help(void);
void merge_files(std::string ifilename, std::string ofilename);
void consume(TaskQueue <int> & q, Options options, PeakIntervals* pintervals, int thread_index);
void fill_queue(const int numcopies, TaskQueue<int> & q);

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
    } else if (PARAMETER_CHECK("-b", 2, parameterLength)) {
      if ((i+1) < argc) {
	options.chipbam = argv[i+1];
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
    } else if (PARAMETER_CHECK("--thread", 8, parameterLength)){
      if ((i+1) < argc) {
    options.n_threads = std::atoi(argv[i+1]);
    i++;
      }
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
  if (options.peakfiletype.empty()) {
    cerr << "****** ERROR: Must specify peakfiletype with -t ******" << endl;
    showHelp = true;
  }

  if (!showHelp) {
    TaskQueue<int> task_queue;
    fill_queue(options.numcopies, task_queue);

    /***************** Main implementation ***************/
    // Perform in bins so we don't keep everything in memory at once
    PeakIntervals* pintervals = \
               new PeakIntervals(options.peaksbed, options.peakfiletype, options.chipbam, options.countindex);

    //BinQueue<GenomeBin> bq;
    //fill_queue(bingenerator, bq);

    // Keeps track of number of reads for naming
    std::vector<std::thread> consumers;
    for (int thread_index=0; thread_index<options.n_threads; thread_index++){
      std::thread cnsmr(std::bind(consume, std::ref(task_queue), options, pintervals, thread_index));
      consumers.push_back(std::move(cnsmr));
    }

    // wait until all threads finish
    for (auto & cnsmr: consumers){
      cnsmr.join();
    }
    
    for (int thread_index=0; thread_index<options.n_threads; thread_index++){
      if (options.paired){
        std::string ifilename_1 = options.outprefix+"/reads_"+std::to_string(thread_index)+"_1.fastq";
        std::string ofilename_1 = options.outprefix+"/reads_1.fastq";
        merge_files(ifilename_1, ofilename_1);

        std::string ifilename_2 = options.outprefix+"/reads_"+std::to_string(thread_index)+"_2.fastq";
        std::string ofilename_2 = options.outprefix+"/reads_2.fastq";
        merge_files(ifilename_2, ofilename_2);
      }else{
        std::string ifilename = options.outprefix+"/reads_"+std::to_string(thread_index)+".fastq";
        std::string ofilename = options.outprefix+"/reads.fastq";
        merge_files(ifilename, ofilename);
      }
    }

    delete pintervals;
    return 0;
    /******************************************************/
  } else {
    simulate_reads_help();
    return 1;
  }
}

/*
 * A thread that keep generate bins
 * */
/*
void fill_queue(BinGenerator bingenerator, BinQueue <GenomeBin> & q){
  while(bingenerator.GotoNextBin()) {
    q.push(bingenerator.GetCurrentBin());
  }
}
*/
void fill_queue(const int numcopies, TaskQueue<int> & q){
  for (int copy_index=0; copy_index<numcopies; copy_index++){
    q.push(copy_index);
  }
}


/*
 * A tread that operate on bins
 * */
void consume(TaskQueue <int> & q, Options options, PeakIntervals* pintervals, int thread_index){
  while (true){
    int copy_index = -1;
    try{
      copy_index = q.pop();
    } catch (std::out_of_range e){
      cerr << e.what() << endl;
      break;
    }

    int total_reads = 0;
    int peakIndexStart = 0;
    int start_offset = 0;
    std::string prev_chrom = "";
    BinGenerator bingenerator(options);
    while (bingenerator.GotoNextBin()){
      // set up
      vector <Fragment> pulldown_fragments, lib_fragments;

      /*** Step 1/2: Shearing + Pulldown ***/
      Pulldown pulldown(options, bingenerator.GetCurrentBin(),\
                            prev_chrom, peakIndexStart, start_offset);
      pulldown.Perform(&pulldown_fragments, pintervals);

      /*** Step 3: Library construction ***/
      LibraryConstructor lc(options);
      lc.Perform(pulldown_fragments, &lib_fragments);

      /*** Step 4: Sequencing ***/
      Sequencer seq(options);
      seq.Sequence(lib_fragments, total_reads, thread_index, copy_index);
    }
  }
}

void merge_files(std::string ifilename, std::string ofilename){
  std::ifstream ifile(ifilename.c_str());
  std::ofstream ofile(ofilename.c_str());
  std::string line;
  while (std::getline(ifile, line)){
    ofile<<line<<"\n";
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
