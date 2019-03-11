#include <cstring>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <sstream>
#include <cmath>
#include <numeric>
#include <time.h>
#include <algorithm>
#include <chrono>

#include "src/bam_io.h"
#include "src/common.h"
#include "src/model.h"
#include "src/options.h"
#include "src/peak_loader.h"
#include "src/fragment.h"
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
bool learn_pcr(const std::string& bamfile, float* geo_rate);
bool learn_frag_paired(const std::string& bamfile, float* alpha, float* beta, bool skip_frag);
bool learn_frag_single(const std::string& bamfile, const std::string& peakfile, const std::string peakfileType,
    const std::int32_t count_colidx, const int intensity_threshold, const int estimate_frag_length,
    float* alpha, float* beta, bool skip_frag);
void search(float& low, float& high, const float mu, const float gs,
    const std::int32_t start_lower_bound, const std::int32_t start_upper_bound, const std::vector<float> start_cdf, const std::vector<float> start_edf,
    const std::int32_t end_lower_bound, const std::int32_t end_upper_bound, const std::vector<float> end_cdf, const std::vector<float> end_edf);
bool calculate_distribution(const std::vector<float> points, const std::int32_t lower_bound, const std::int32_t upper_bound,
    std::vector<float>& pdf, std::vector<float>& cdf, std::vector<float>& edf);
float calculate_gs(const std::int32_t lower_bound, const std::int32_t upper_bound,
    const std::vector<float> cdf, const std::vector<float> edf, const std::vector<float> pdf, const float mu,
    const float offset, const std::int32_t bound_baseline);
float calculate_eexpl(const float k, const float theta, const std::int32_t lower_bound, const std::int32_t upper_bound,
    const std::vector<float> cdf, const std::vector<float> edf, const float mu);
float calculate_edf(const float x, const std::int32_t lower_bound, const std::int32_t upper_bound,
    const std::vector<float> cdf, const std::vector<float> edf, const float mu);
float calculate_gamma_pdf(const float x, const float k, const float theta);
bool learn_frag_paired(const std::string& bamfile, float* alpha, float* beta);


bool learn_frag_paired(const std::string& bamfile, float* alpha, float* beta, bool skip_frag) {
  /*
    Learn fragment length distribution from an input BAM file
    Fragment lengths follow a gamma distribution.

    Inputs:
    - bamfile (std::string): path to the BAM file
    - skip_frag (bool): Skip fragment parameter learning process
    Outputs:
    - alpha (float): parameter of gamma distribution
    - beta  (float): parameter of gamma distribution
   */

  /* Learning process should be skipped */
  if (skip_frag){return true;}

  /* First, get a vector of the fragment lengths */
  int maxreads = 10000; int numreads = 0; // don't look at more than this many reads
  BamCramReader bamreader(bamfile);
  const BamHeader* bamheader = bamreader.bam_header();
  // Get first chrom to look at fragment lengths
  std::vector<std::string> seq_names = bamheader->seq_names();
  std::vector<std::uint32_t> seq_lengths = bamheader->seq_lengths();
  std::vector<std::int32_t> fraglengths;
  std::int32_t tlen;
  BamAlignment aln;

  // Randomly grab fragments
  if (seq_names.size() > 0 && seq_lengths.size() > 0) {
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    srand(seed);

    int guard_count = 0; // prevent the dead loop
    while(numreads<maxreads && guard_count<(100*maxreads)) {
      int chrom = rand() % seq_names.size();
      int start = rand() % seq_lengths[chrom];
      
      if ((seq_names[chrom].find("_") != std::string::npos) || (seq_names[chrom] == "chrM")){continue;}
      bamreader.SetRegion(seq_names[chrom], start, seq_lengths[chrom]);

      if (!bamreader.GetNextAlignment(aln)) {continue;}
      if (aln.IsDuplicate()) {continue;}
      if ( (!aln.IsMapped()) || aln.IsSecondary()){continue;}
      tlen = aln.TemplateLength();
      if (tlen > 0) {
        fraglengths.push_back(abs(tlen));
        numreads++;
      }
      guard_count++;
      //    cerr << abs(tlen) << endl; // if you want to print out for debugging
    }

    // find median to use to filter
    nth_element(fraglengths.begin(), fraglengths.begin() + fraglengths.size()/2, fraglengths.end());
    int first = fraglengths[fraglengths.size()/2];
    nth_element(fraglengths.begin(), fraglengths.begin() + fraglengths.size()/2 - 1, fraglengths.end());
    int second = fraglengths[fraglengths.size()/2 - 1];
    float median = (first + second)/2;

    for (int fragment=0; fragment < fraglengths.size(); fragment++)
    {
      if (fraglengths[fragment] > (3*median))
      {
        fraglengths.erase(fraglengths.begin() + fragment);
        fragment--;
      }
    }
  } else {
    return false;
  }

  float total_frag_len = 0;     // sum of all the frag lengths

  // get sum of all frag lengths and log sum of each frag length
  for (int frag = 0; frag < fraglengths.size(); frag++)
  {
    total_frag_len += fraglengths[frag];
  }

  // mean of frag lengths and log mean of fraglengths
  float mean_frag_length = total_frag_len/fraglengths.size();

  /* Using method of moments to estimate the shape and scale parameters */
  float moment_sum = 0;
  for (int n = 0; n < fraglengths.size(); n++)
  {
    moment_sum += (fraglengths[n] - mean_frag_length)*(fraglengths[n] - mean_frag_length);
  }
  *beta = (moment_sum/(mean_frag_length*fraglengths.size()));
  *alpha = mean_frag_length / *beta;

  if (DEBUG) { 
    std::stringstream ss;
    ss << "Learned fragment length params alpha: " << *alpha << " and beta: " << *beta;
    PrintMessageDieOnError(ss.str(), M_DEBUG);
  }
  return true;
}


bool learn_frag_single(const std::string& bamfile,
                        const std::string& peakfile, const std::string peakfileType,
                        const std::int32_t count_colidx, const int intensity_threshold,
                        const int estimate_frag_length,
                        float* alpha, float* beta, bool skip_frag) {
  /*
    Predict fragment length distribution from an input BAM file (single-end reads)
    Fragment lengths follow a gamma distribution.

    Inputs:
    - bamfile (std::string): path to the BAM file
    Outputs:
    - alpha (float): parameter of gamma distribution
    - beta  (float): parameter of gamma distribution
   */

  /* Learning process should be skipped */
  if (skip_frag){return true;}

  /* First read peaks and restrict reads origins */
  // peak intensity is only used for filtering out unreliable reads
  PeakLoader peakloader(peakfile, peakfileType, "", count_colidx);
  std::vector<Fragment> peaks;
  if (!peakloader.Load(peaks)) PrintMessageDieOnError("Error loading peaks from " + peakfile, M_ERROR);
  for (int peak_idx=peaks.size()-1;peak_idx>=0;peak_idx--){
    if (peaks[peak_idx].score < intensity_threshold) peaks.erase(peaks.begin()+peak_idx);
  }

  /* Read reads from the BAM file */
  BamCramReader bamreader(bamfile);
  std::vector<float> starts;
  std::vector<float> ends;
  for (int peak_index=0; peak_index<peaks.size(); peak_index++){
    bamreader.SetRegion(peaks[peak_index].chrom, peaks[peak_index].start, peaks[peak_index].start+peaks[peak_index].length);
    BamAlignment aln;
    std::vector<float> starts_in_peak;
    std::vector<float> ends_in_peak;
    while (bamreader.GetNextAlignment(aln)){
      if (aln.IsDuplicate()) {continue;}
      if ( (!aln.IsMapped()) || aln.IsSecondary()){continue;}
      float aln_start = aln.Position();
      float aln_end = aln.GetEndPosition();

      if(aln.IsReverseStrand()){
        if ( ((aln_end-estimate_frag_length) >= peaks[peak_index].start)
                &&(aln_end <= peaks[peak_index].start+peaks[peak_index].length)){
          ends_in_peak.push_back(aln_end);
        }
      }else{
        if ((aln_start >= peaks[peak_index].start) 
                &&( (aln_start+estimate_frag_length) <= peaks[peak_index].start+peaks[peak_index].length)){
          starts_in_peak.push_back(aln_start);
        }
      }
    }
    
    if ((starts_in_peak.size() == 0) || (ends_in_peak.size() == 0)) continue;
    float avg_start = std::accumulate(starts_in_peak.begin(), starts_in_peak.end(), 0.0)/ (float) starts_in_peak.size();
    float avg_end = std::accumulate(ends_in_peak.begin(), ends_in_peak.end(), 0.0)/ (float) ends_in_peak.size();
    float mid_pos = (avg_start + avg_end) / 2.0;

    // note:
    // starts will be float variables due to "mid_pos"
    // starts are less then zero in most of the time; ends are the otherwise
    for (int start_idx=0; start_idx<starts_in_peak.size(); start_idx++) starts.push_back(starts_in_peak[start_idx] - mid_pos);
    for (int end_idx=0; end_idx<ends_in_peak.size(); end_idx++) ends.push_back(ends_in_peak[end_idx] - mid_pos);
  }

  /* mean value of fragment length*/
  float mean_frag_length = std::accumulate(ends.begin(), ends.end(), 0.0)/ (float) ends.size()
                        - std::accumulate(starts.begin(), starts.end(), 0.0)/ (float) starts.size();

  /* calculate CDF of starts and ends */
  /* start */
  std::int32_t start_lower_bound = std::floor(*std::min_element(starts.begin(), starts.end()));
  std::int32_t start_upper_bound = std::ceil(*std::max_element(starts.begin(), starts.end()));
  vector<float> start_pdf(start_upper_bound-start_lower_bound+1, 0);
  vector<float> start_cdf(start_upper_bound-start_lower_bound+1, 0);
  vector<float> start_edf(start_upper_bound-start_lower_bound+1, 0);
  if (!calculate_distribution(starts, start_lower_bound, start_upper_bound,
                        start_pdf, start_cdf, start_edf)){
    PrintMessageDieOnError("Error happened when calculating position distributions", M_ERROR);
  }

  std::int32_t end_lower_bound = std::floor(*std::min_element(ends.begin(), ends.end()));
  std::int32_t end_upper_bound = std::ceil(*std::max_element(ends.begin(), ends.end()));
  vector<float> end_pdf(end_upper_bound-end_lower_bound+1, 0);
  vector<float> end_cdf(end_upper_bound-end_lower_bound+1, 0);
  vector<float> end_edf(end_upper_bound-end_lower_bound+1, 0);
  if (!calculate_distribution(ends, end_lower_bound, end_upper_bound,
                        end_pdf, end_cdf, end_edf)){
    PrintMessageDieOnError("Error happened when calculating position distributions", M_ERROR);
  }

  float gs1 = calculate_gs(end_lower_bound, end_upper_bound, end_cdf, end_edf,
                            start_pdf, mean_frag_length,
                            (mean_frag_length + mean_frag_length/3.0)/2.0, start_lower_bound);
  float gs2 = calculate_gs(start_lower_bound, start_upper_bound, start_cdf, start_edf,
                            end_pdf, -mean_frag_length,
                            -(mean_frag_length + mean_frag_length/3.0)/2.0, end_lower_bound);
  float gs = std::max(gs1, gs2);
  float low = 200;
  float high = 8000;
  for(int experiment_idx=0; experiment_idx<35; experiment_idx++){
    search(low, high, mean_frag_length, gs,
            start_lower_bound, start_upper_bound, start_cdf, start_edf,
            end_lower_bound, end_upper_bound, end_cdf, end_edf);
  }
  float var_frag_length = (low + high) / 2.0;

  // mean_frag_length, var_frag_length
  // --> k and theta
  *beta = var_frag_length / mean_frag_length;
  *alpha = mean_frag_length / *beta;
  if (DEBUG) {
    std::stringstream ss;
    ss << "Learned fragment length params alpha: " << *alpha << " and beta: " << *beta;
    PrintMessageDieOnError(ss.str(), M_DEBUG);
  }
  return true;
}

void search(float& low, float& high, const float mu, const float gs,
            const std::int32_t start_lower_bound, const std::int32_t start_upper_bound,
            const std::vector<float> start_cdf, const std::vector<float> start_edf,
            const std::int32_t end_lower_bound, const std::int32_t end_upper_bound,
            const std::vector<float> end_cdf, const std::vector<float> end_edf){
  // calculate the score for the current setting
  float shape = mu * mu / (float(low + high)/2.0);
  float scale = (float(low + high)/2.0) / mu;
  float eexpl1 = calculate_eexpl(shape, scale, start_lower_bound, start_upper_bound, start_cdf, start_edf, -mu);
  float eexpl2 = calculate_eexpl(shape, scale, end_lower_bound, end_upper_bound, end_cdf, end_edf, mu);
  float eexpl = std::min(eexpl1, eexpl2);

  if ((eexpl-gs) > 0){
    high = (float(low + high)/2.0);
  }else{
    low = (float(low + high)/2.0);
  }
  return;
}


bool calculate_distribution(const std::vector<float> points,
        const std::int32_t lower_bound, const std::int32_t upper_bound,
        std::vector<float>& pdf, std::vector<float>& cdf, std::vector<float>& edf){
  // pdf
  for(int pt_idx=0; pt_idx<points.size(); pt_idx++){
    std::uint32_t key = std::ceil(points[pt_idx]) - lower_bound;
    pdf[key] += 1;
  }
  for(int pdf_idx=0; pdf_idx<pdf.size(); pdf_idx++) pdf[pdf_idx] /= (float) points.size(); // norm

  // cdf
  cdf[0] = pdf[0];
  for(int cdf_idx=1; cdf_idx<cdf.size(); cdf_idx++){
    cdf[cdf_idx] = cdf[cdf_idx-1] + pdf[cdf_idx];
  }

  // edf - left
  std::vector<float> edf_l(upper_bound-lower_bound+1, 0);
  for(int pt_idx=0; pt_idx<points.size(); pt_idx++){
    std::uint32_t key = std::floor(points[pt_idx]) - lower_bound;
    edf_l[key] += (points[pt_idx] - std::floor(points[pt_idx]));
  }
  for(int edf_idx=0; edf_idx<edf_l.size(); edf_idx++){
    edf_l[edf_idx] /= (float) points.size(); // norm
  }
  for(int edf_idx=edf_l.size()-1; edf_idx>0; edf_idx--){
    edf_l[edf_idx-1] += edf_l[edf_idx] + (1 - cdf[edf_idx]);
  }

  // edf - right
  vector<float> edf_r(upper_bound-lower_bound+1, 0);
  for(int pt_idx=0; pt_idx<points.size(); pt_idx++){
    std::uint32_t key = std::ceil(points[pt_idx]) - lower_bound;
    edf_r[key] += (std::ceil(points[pt_idx]) - points[pt_idx]);
  }
  for(int edf_idx=0; edf_idx<edf_r.size(); edf_idx++){
    edf_r[edf_idx] /= (float) points.size(); // norm
  }
  for(int edf_idx=0; edf_idx<(edf_r.size()-1); edf_idx++){
    edf_r[edf_idx+1] += (edf_r[edf_idx] + cdf[edf_idx]);
  }
  // edf-normed
  for (int edf_idx=0; edf_idx<edf.size(); edf_idx++){
    edf[edf_idx] = (edf_l[edf_idx]+edf_r[edf_idx]);
  }

  return true;
}

float calculate_gs(const std::int32_t lower_bound, const std::int32_t upper_bound,
        const std::vector<float> cdf, const std::vector<float> edf, 
        const std::vector<float> pdf, const float mu,
        const float offset, const std::int32_t bound_baseline){
  float gs = 0;
  for(int pdf_idx=0; pdf_idx<pdf.size(); pdf_idx++){
    int abs_pos = pdf_idx + bound_baseline; // relative_pos + lower_bound - the_other_lower_bound
    gs += (pdf[pdf_idx] * calculate_edf((abs_pos/3.0 + offset), lower_bound, upper_bound, cdf, edf, mu));
  }
  return gs;
}

float calculate_eexpl(const float k, const float theta,
                        const std::int32_t lower_bound, const std::int32_t upper_bound,
                        const std::vector<float> cdf, const std::vector<float> edf, const float mu){
  float result =0;
  for(int idx=0; idx<800; idx++){
    float edf_idx = idx;
    if (mu < 0){
      edf_idx = -edf_idx;
    }
    result += (calculate_gamma_pdf(idx, k, theta) *
            calculate_edf(edf_idx/2.0, lower_bound, upper_bound, cdf, edf, mu));
  }
  return result;
}

float calculate_edf(const float x,
                    const std::int32_t lower_bound, const std::int32_t upper_bound,
                    const std::vector<float> cdf, const std::vector<float> edf, const float mu){
  if (x <= lower_bound){
    return (mu/2.0 - x);
  } else if (x > upper_bound){
    return (x - mu/2.0);
  } else {
    float key = std::floor(x-lower_bound);
    return (edf[key] + (x-std::floor(x)) * (2 * cdf[key] - 1));
  }
}

float calculate_gamma_pdf(const float x, const float k, const float theta){
  if (x <= 0 || k<= 0 || theta<=0){
    return 0.0;
  }else{
    return (exp(-x/theta) * pow(x, k - 1.0) / (pow(theta, k) * std::tgamma(k)));
  }
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
  for(int peak_index = 0; peak_index < peaks.size(); peak_index++){
    plen += peaks[peak_index].length;
  }

  // calculate f and s, and then ab_ratio
  float f = (float)plen / (float) (peakloader.total_genome_length);
  float s = (float)(peakloader.tagcount_in_peaks) / (float)(peakloader.total_tagcount);

  // outputs
  float ab_ratio = (s/(1-s)) * ((1-f)/f);
  *ab_ratio_ptr = ab_ratio;
  *f_ptr = f;
  *s_ptr = s;
  if (DEBUG) {
    std::stringstream ss;
    ss <<"peak-length: "<<plen<<"\ttotal: "<< peakloader.total_genome_length<<"\n";
    ss <<"#reads in peaks: "<<peakloader.tagcount_in_peaks<<"\t#reads: "<<peakloader.total_tagcount<<"\n";
    ss <<"ab_ratio: " << *ab_ratio_ptr << " f: " << *f_ptr << " s: " << *s_ptr << "\n";
    PrintMessageDieOnError(ss.str(), M_DEBUG);
  }
  return true;
}

bool learn_pcr(const std::string& bamfile, float* geo_rate){
  BamCramReader bamreader(bamfile);
  const BamHeader* bamheader = bamreader.bam_header();
  std::vector<std::string> seq_names = bamheader->seq_names();
  std::vector<uint32_t> seq_lengths = bamheader->seq_lengths();

  std::map<int, int> n_pcr_copies = {{1, 0}};
  std::map<std::string, int> duplicated_reads;
  int count_unmapped = 0;
  for (int seq_index=0; seq_index<seq_names.size(); seq_index++){
    if ((seq_names[seq_index].find("_") == std::string::npos) && (seq_names[seq_index] != "chrM")){
      bamreader.SetRegion(seq_names[seq_index], 0, seq_lengths[seq_index]);
      BamAlignment aln;
      while (bamreader.GetNextAlignment(aln)){
       if ( (!aln.IsMapped()) || aln.IsSecondary()){count_unmapped+=1; continue;}
       if (aln.IsDuplicate()){
        std::string read_location = seq_names[seq_index] + ":" + std::to_string(aln.Position());
        if (duplicated_reads.find(read_location) != duplicated_reads.end()){
          duplicated_reads[read_location] += 1;
        }else{
          duplicated_reads[read_location] = 2;
        }
       }else{
        n_pcr_copies[1] += 1;
       }
      }
    }
  }

  //for (std::map<Fragment, int, bool(*)(Fragment, Fragment)>::iterator it=duplicated_reads.begin();
  for (std::map<std::string, int>::iterator it=duplicated_reads.begin();
          it!=duplicated_reads.end(); ++it){
    if (n_pcr_copies.find(it->second) != n_pcr_copies.end()){
      n_pcr_copies[it->second] += 1;
    }else{
      n_pcr_copies[it->second] = 1;
    }
  }

  int total_n_samples = 0;
  *geo_rate = 0;
  for (std::map<int, int>::iterator it=n_pcr_copies.begin(); it!=n_pcr_copies.end(); ++it){
    *geo_rate += (it->first) * (it->second);
    total_n_samples += it->second;
  }
  *geo_rate /= float(total_n_samples);
  *geo_rate = 1.0/(*geo_rate);
  if (DEBUG) {
    std::stringstream ss;
    ss << "Geo_rate: " << *geo_rate;
    PrintMessageDieOnError(ss.str(), M_DEBUG);
  }
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
    	options.remove_pct = std::atof(argv[i+1]);
    	i++;
      }
    } else if (PARAMETER_CHECK("--skip-frag", 11, parameterLength)){
      if ((i) < argc) {
    options.skip_frag = true;
    i++;
      }
    } else if (PARAMETER_CHECK("--thres", 7, parameterLength)){
      if ((i+1) < argc){
    options.intensity_threshold = std::atoi(argv[i+1]);
    i++;
      }
    } else if (PARAMETER_CHECK("--paired", 8, parameterLength)){
      options.paired = true;
    } else if (PARAMETER_CHECK("--est", 5, parameterLength)){
      if ((i+1) < argc){
    options.estimate_frag_length = std::atoi(argv[i+1]);
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
    ChIPModel model;

    /*** Learn fragment size disbribution parameters ***/
    float frag_param_a = -1;
    float frag_param_b = -1;
    if (options.paired){
      if (!learn_frag_paired(options.chipbam, &frag_param_a, &frag_param_b, options.skip_frag)) {
        PrintMessageDieOnError("Error learning fragment length distribution", M_ERROR);
      }
    }else{
      if (!learn_frag_single(options.chipbam, options.peaksbed, options.peakfiletype,
                  options.countindex, options.intensity_threshold, options.estimate_frag_length,
                  &frag_param_a, &frag_param_b, options.skip_frag)) {
        PrintMessageDieOnError("Error learning fragment length distribution", M_ERROR);
      }
    }
    model.SetFrag(frag_param_a, frag_param_b);

    /*** Learn pulldown ratio parameters ***/
    float ab_ratio;
    float s = -1;
    float f = -1;
    if (!learn_ratio(options.chipbam, options.peaksbed, options.peakfiletype, options.countindex, options.remove_pct,
        &ab_ratio, &s, &f)){
      PrintMessageDieOnError("Error learning pulldown ratio", M_ERROR);
    }
    model.SetF(f);
    model.SetS(s);

    /*** Learn PCR geometric distribution parameter **/
    float geo_rate = -1;
    if (!learn_pcr(options.chipbam, &geo_rate)){
      PrintMessageDieOnError("Error learning PCR rate", M_ERROR);
    }
    model.SetPCR(geo_rate);

    /*** Output learned model **/
    model.WriteModel(options.outprefix + ".json");
    model.PrintModel();
    return 0;
  } else {
    learn_help();
    return 1;
  }
}


void learn_help(void) {
  Options options;
  cerr << "\nTool:    chipmunk learn" << endl;
  cerr << "Version: " << _GIT_VERSION << "\n";    
  cerr << "Summary: Learn parameters from a ChIP dataset." << endl << endl;
  cerr << "Usage:   " << PROGRAM_NAME << " learn -b reads.bam -p peak.bed -o outprefix [OPTIONS] " << endl << endl;
  cerr << "[Required arguments]: " << "\n";
  cerr << "         -b <reads.bam>:     BAM file with ChIP reads (.bai index required)" << "\n";
  cerr << "         -p <peaks.bed>:     BED file with peak regions (Homer format or BED format)" << "\n";
  cerr << "         -t <peakfile_type>: File type of the input peak file. Only `homer` or `bed` supported." << "\n";
  cerr << "         -o <outprefix>:     Prefix for output files" << "\n";
  cerr << "         -c <int>:           The index of the BED file column used to score each peak (index starting from 1)" << "\n";
  cerr << "[Optional arguments]: " << "\n";
  cerr << "         -r <float>:         Ratio of high score peaks to ignore\n"
       << "                             Default: " <<options.remove_pct<< "\n";
  cerr << "         --est <int>:        Estimated fragment length\n"
       << "                             Default: " <<options.estimate_frag_length<< "\n";
  cerr << "[BAM-file arguments]: " << "\n";
  cerr << "         --paired:           Loading paired-end reads\n"
       << "                             Default: false\n";
  cerr << "         --thres <float>:    Threshold for peak scores in single-end read length estimation\n"
       << "                             Default: " <<options.intensity_threshold<< "\n";
  cerr << "\n";
  cerr  << "[ General help ]:" << endl;
  cerr  << "    --help        "  << "Print this help menu.\n";
  cerr << "\n";
  exit(1);
}
