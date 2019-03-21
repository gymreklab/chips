#include <cstring>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <stdio.h>
#include <random>

#include "src/bingenerator.h"
#include "src/common.h"
#include "src/fragment.h"
#include "src/library_constructor.h"
#include "src/model.h"
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
void consume(TaskQueue <int> & q, Options options, PeakIntervals* pintervals, const std::vector<int>& reads_per_copy, int thread_index);
void fill_queue(const int numcopies, TaskQueue<int> & q);
void GetReadsPerCopy(std::vector<int>* reads_per_copy, const Options& options);

int simulate_reads_main(int argc, char* argv[]) {
  bool showHelp = false;
  Options options;

  // keep track of model params. Command-line options override json
  ChIPModel model;

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
    } else if (PARAMETER_CHECK("-v", 2, parameterLength)) {
      options.verbose = true;
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
	options.gamma_k = atof(params[0].c_str());
	options.gamma_theta = atof(params[1].c_str());
	model.SetFrag(options.gamma_k, options.gamma_theta);
	i++;
      }
    } else if (PARAMETER_CHECK("--spot", 6, parameterLength)) {
      if ((i+1) < argc) {
	options.ratio_s = atof(argv[i+1]);
	model.SetS(options.ratio_s);
	i++;
      }
    } else if (PARAMETER_CHECK("--frac", 6, parameterLength)) {
      if ((i+1) < argc) {
	options.ratio_f = atof(argv[i+1]);
	model.SetF(options.ratio_f);
	i++;
      }
    } else if (PARAMETER_CHECK("--recomputeF", 12, parameterLength)) {
      options.recompute_f = true;
    } else if (PARAMETER_CHECK("--model", 7, parameterLength)) {
      if ((i+1) < argc) {
	options.model_file = argv[i+1];
	model.ReadFromJSON(options.model_file);
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
    } else if (PARAMETER_CHECK("--noscale", 9, parameterLength)) {
      options.noscale = true;
    } else if (PARAMETER_CHECK("--thread", 8, parameterLength)){
      if ((i+1) < argc) {
	options.n_threads = std::atoi(argv[i+1]);
	i++;
      }
    } else if (PARAMETER_CHECK("--sequencer", 11, parameterLength)){
      if ((i+1) < argc){
	options.sequencer_type = argv[i+1];
	i++;
      }
    } else if (PARAMETER_CHECK("--sub", 5, parameterLength)){
      if ((i+1) < argc){
	options.sub_rate = std::atof(argv[i+1]);
	i++;
      }
    } else if (PARAMETER_CHECK("--ins", 5, parameterLength)){
      if ((i+1) < argc){
	options.ins_rate = std::atof(argv[i+1]);
	i++;
      }
    } else if (PARAMETER_CHECK("--del", 5, parameterLength)){
      if ((i+1) < argc){
	options.del_rate = std::atof(argv[i+1]);
	i++;
      }
    } else if (PARAMETER_CHECK("--pcr_rate", 10, parameterLength)){
      if ((i+1) < argc){
	options.pcr_rate = std::atof(argv[i+1]);
	model.SetPCR(options.pcr_rate);
	i++;
      }
    } else {
      cerr << endl << "*****ERROR: Unrecognized parameter: " << argv[i] << " *****" << endl << endl;
      showHelp = true;
    }
  }
  // Set model defaults if not set
  model.UpdateParams(options);
  // Keep final model in options since other classes use that
  model.UpdateOptions(options);

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
    // Print out parsed model
    PrintMessageDieOnError("Running simulate with the following model", M_PROGRESS);
    model.PrintModel();

    // Set up jobs
    TaskQueue<int> task_queue;
    fill_queue(options.numcopies, task_queue);

    // Determine number of reads per copy
    std::vector<int> reads_per_copy;
    GetReadsPerCopy(&reads_per_copy, options);

    // Remove previous existing fastqs
    if (options.paired){
      // read - first pair
      std::string reads1 = options.outprefix+"_1.fastq";
      std::remove(reads1.c_str());
      for (int thread_index=0; thread_index<options.n_threads; thread_index++){
        reads1 = options.outprefix+"_"+std::to_string(thread_index)+"_1.fastq";
        std::remove(reads1.c_str());
      }
      // read - second pair
      std::string reads2 = options.outprefix+"_2.fastq";
      std::remove(reads2.c_str());
      for (int thread_index=0; thread_index<options.n_threads; thread_index++){
        reads2 = options.outprefix+"_"+std::to_string(thread_index)+"_2.fastq";
        std::remove(reads2.c_str());
      }
    }
    else{
      std::string reads = options.outprefix+".fastq";
      std::remove(reads.c_str());
      for (int thread_index=0; thread_index<options.n_threads; thread_index++){
        reads = options.outprefix+"_"+std::to_string(thread_index)+".fastq";
        std::remove(reads.c_str());
      }
    }

    /***************** Main implementation ***************/
    // Perform in bins so we don't keep everything in memory at once
    PrintMessageDieOnError("Loading the input ChIP-seq peak file (and BAM file if given)", M_PROGRESS);
    PeakIntervals* pintervals = \
               new PeakIntervals(options, options.peaksbed, options.peakfiletype, options.chipbam, options.countindex);
    if (options.recompute_f) {
      RefGenome ref_genome(options.reffa);
      float f = pintervals->total_bound_length/ref_genome.GetGenomeLength(); // TODO recompute with new frag score
      if (f<0) {
	PrintMessageDieOnError("Error. Estimated --frac negative. Likely overflow error on genome size", M_ERROR);
      }
      model.SetF(f);
      model.UpdateOptions(options);
      PrintMessageDieOnError("Recomputed --frac. New model:", M_PROGRESS);
      model.PrintModel();
    }

    // Create threads
    PrintMessageDieOnError("Simulating reads based on the input profile", M_PROGRESS);
    std::vector<std::thread> consumers;
    for (int thread_index=0; thread_index<options.n_threads; thread_index++){
      std::thread cnsmr(std::bind(consume, std::ref(task_queue), options, pintervals, reads_per_copy, thread_index));
      consumers.push_back(std::move(cnsmr));
    }

    // wait until all threads finish
    for (auto & cnsmr: consumers){
      cnsmr.join();
    }
   
 
    PrintMessageDieOnError("Writing reads into file", M_PROGRESS);
    for (int thread_index=0; thread_index<options.n_threads; thread_index++){
      if (options.paired){
        std::string ifilename_1 = options.outprefix+"_"+std::to_string(thread_index)+"_1.fastq";
        std::string ofilename_1 = options.outprefix+"_1.fastq";
        merge_files(ifilename_1, ofilename_1);
        std::remove(ifilename_1.c_str());

        std::string ifilename_2 = options.outprefix+"_"+std::to_string(thread_index)+"_2.fastq";
        std::string ofilename_2 = options.outprefix+"_2.fastq";
        merge_files(ifilename_2, ofilename_2);
        std::remove(ifilename_2.c_str());
      }else{
        std::string ifilename = options.outprefix+"_"+std::to_string(thread_index)+".fastq";
        std::string ofilename = options.outprefix+".fastq";
        merge_files(ifilename, ofilename);
        std::remove(ifilename.c_str());
      }
    }

    delete pintervals;
    PrintMessageDieOnError("Done!", M_PROGRESS);
    return 0;
    /******************************************************/
  } else {
    simulate_reads_help();
    return 1;
  }
}

/*
 * Determine the number of reads per genome copy
 * */
void GetReadsPerCopy(std::vector<int>* reads_per_copy, const Options& options) {
  // Initialize vector
  reads_per_copy->clear();
  for (size_t i=0; i<options.numcopies; i++) {
    reads_per_copy->push_back(0);
  }
  // Assign each read to a copy
  std::random_device rd;
  std::mt19937 rng(rd());
  std::uniform_int_distribution<int> uni(0, options.numcopies-1);
  int copynum;
  for (size_t i=0; i<options.numreads; i++) {
    copynum = uni(rng);
    (*reads_per_copy)[copynum] += 1;
  }
}

/*
 * A thread that generates genome copies
 * */
void fill_queue(const int numcopies, TaskQueue<int> & q){
  for (int copy_index=0; copy_index<numcopies; copy_index++){
    q.push(copy_index);
  }
}


/*
 * A thread that operate on a single genome copy
 * */
void consume(TaskQueue <int> & q, Options options, PeakIntervals* pintervals, const std::vector<int>& reads_per_copy, int thread_index){
  while (true){
    int copy_index = -1;
    try{
      copy_index = q.pop();
    } catch (std::out_of_range e){
      //cerr << e.what() << endl;
      break;
    }

    if ((copy_index > 0) && (copy_index%100 == 0)) {
        int job_percentage = (int) (100 * copy_index / (float) options.numcopies);
        PrintMessageDieOnError("Simulated " + std::to_string(job_percentage) +"% reads.", M_PROGRESS);
    }

    if (reads_per_copy[copy_index] == 0) {
      continue; // If we're not going to get any reads, don't bother simulating
    }
    int total_reads = 0;
    int peakIndexStart = 0;
    int start_offset = 0;
    std::string prev_chrom = "";
    BinGenerator bingenerator(options);
    // set up. Clear pulldown each time. Append to lib_fragments and sequence at the end
    vector <Fragment> pulldown_fragments, lib_fragments;
    while (bingenerator.GotoNextBin()){
      if (options.verbose) {
	stringstream ss;
	ss << "Processing bin " << bingenerator.GetCurrentBinStr() << " " << copy_index;
	PrintMessageDieOnError(ss.str(), M_PROGRESS);
      }

      /*** Step 1/2: Shearing + Pulldown ***/
      Pulldown pulldown(options, bingenerator.GetCurrentBin(),\
			prev_chrom, peakIndexStart, start_offset);
      pulldown.Perform(&pulldown_fragments, pintervals);

      /*** Step 3: Library construction ***/
      LibraryConstructor lc(options);
      lc.Perform(pulldown_fragments, &lib_fragments);

      /*** Cleanup for next bin ***/
      pulldown_fragments.clear();
    }
    /*** Step 4: Sequencing ***/
    Sequencer seq(options);
    seq.Sequence(lib_fragments, reads_per_copy[copy_index], total_reads, thread_index, copy_index);
  }
}

void merge_files(std::string ifilename, std::string ofilename){
  std::ifstream ifile(ifilename.c_str());
  std::ofstream ofile(ofilename.c_str(), std::ofstream::app);
  std::string line;
  while (std::getline(ifile, line)){
    ofile<<line<<"\n";
  }
}

void simulate_reads_help(void) {
  Options options;
  cerr << "\nTool:    chipmunk simreads" << endl;
  cerr << "Version: " << _GIT_VERSION << "\n";    
  cerr << "Summary: Simulate ChIP-seq reads for a set of peaks." << endl << endl;
  cerr << "Usage:   " << PROGRAM_NAME << " simreads -p peaks.bed -f ref.fa -o outprefix [OPTIONS] " << endl;
  cerr << "\n[Required arguments]: " << "\n";
  cerr << "     -p <peaks.bed>: BED file with peak regions" << "\n";
  cerr << "     -t <str>: The file format of your input peak file. Only `homer` or `bed` are supported." << "\n";
  cerr << "     -f <ref.fa>: FASTA file with reference genome" << "\n";
  cerr << "     -o <outprefix>: Prefix for output files" << "\n";
  cerr << "\n[Experiment parameters]: " << "\n";
  cerr << "     --numcopies <int>: Number of copies of the genome to simulate\n"
       << "                        Default: " << options.numcopies << "\n";
  cerr << "     --numreads <int> : Number of reads (or read pairs) to simulate\n"
       << "                        Default: " << options.numreads << "\n";
  cerr << "     --readlen <int>  : Read length to generate\n"
       << "                        Default: " << options.readlen << "\n";
  cerr << "     --paired         : Simulate paired-end reads\n"
       << "                        Default: false \n";
  cerr << "\n[Model parameters]: " << "\n";
  cerr << "     --model <str>               : JSON file with model parameters (e.g. from running learn\n";
  cerr << "                                   Setting parameters below overrides anything in the JSON file\n";
  cerr << "     --gamma-frag <float>,<float>: Parameters for fragment length distribution (alpha, beta).\n"
       << "                                   Default: " << options.gamma_k << ","
       << options.gamma_theta << "\n";
  cerr << "     --spot <float>              : SPOT score (fraction of reads in peaks) \n"
       << "                                   Default: " << options.ratio_s << "\n";
  cerr << "     --frac <float>              : Fraction of the genome that is bound \n"
       << "                                   Default: " << options.ratio_f << "\n";
  cerr << "     --recomputeF                : Recompute --frac param based on input peaks.\n";
  cerr << "     --pcr_rate <float>          : The rate of geometric distribution for PCR simulation\n"
       << "                                   Default: " << options.pcr_rate << "\n";
  cerr << "\n[Peak scoring: choose one]: " << "\n";
  cerr << "     -b <reads.bam>              : Read BAM file used to score each peak\n"
       << "                                 : Default: None (use the scores from the peak file)\n";
  cerr << "     -c <int>                    : The index of the BED file column used to score each peak (index starting from 1). Required if -b not used.\n"
       << "                                 : Default: " << options.countindex << "\n";
  cerr << "     --noscale                   : Don't scale peak scores by the max score.\n";                   
  cerr << "\n[Other options]: " << "\n";
  cerr << "     --region <str>              : Only simulate reads from this region chrom:start-end\n"
       << "                                   Default: genome-wide \n";
  cerr << "     --binsize <int>             : Consider bins of this size when simulating\n"
       << "                                 : Default: " << options.binsize << "\n";
  cerr << "     --thread <int>              : Number of threads used for computing\n"
       << "                                 : Default: " << options.n_threads << "\n";
  cerr << "     --sequencer <std>           : Sequencing error values\n"
       << "                                 : Default: None (no sequencing errors)\n";
  cerr << "     --sub <float>               : Customized substitution value in sequecing\n";
  cerr << "     --ins <float>               : Customized insertion value in sequecing\n";
  cerr << "     --del <float>               : Customized deletion value in sequecing\n";
  cerr << "\n";
  cerr  << "[ General help ]:" << endl;
  cerr  << "    --help        "  << "Print this help menu.\n";
  exit(1);
}
