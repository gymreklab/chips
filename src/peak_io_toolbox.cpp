// Peak Loader Tool Box
//

#include "src/peak_io_toolbox.h"

PeakReader::PeakReader(const std::string& _peakfile){
  peakfile = _peakfile;
}

bool PeakReader::HomerPeakReader(std::vector<Fragment>& peaks,
        const std::int32_t count_colidx, const std::string region, const bool noscale){
  std::ifstream infile(peakfile.c_str());
  if(infile.fail()) PrintMessageDieOnError("Input path \"" + peakfile + "\" does not exist", M_ERROR);

  std::string line;
  while (std::getline(infile, line)){
    if (line.find('#') != 0){
      std::string chr;
      std::int32_t start;
      std::int32_t end;
      std::int32_t length;
      float count = -1;
      
      std::stringstream linestream(line);
      std::string element;
      uint32_t elem_idx=0;
      while(!linestream.eof()){
        linestream >> element;
        if (elem_idx == 1){
          chr = element;
        }else if(elem_idx == 2){
          start = std::stol(element);
        }else if(elem_idx == 3){
          end = std::stol(element);
        }else if((count_colidx > 0) && (elem_idx == count_colidx-1)){
          count = std::stof(element);
        }
        elem_idx++;
      }
      if (region.empty()){
        length = end-start;
        Fragment peak_location(chr, start, length, count);
        peaks.push_back(peak_location);
      }else{
        std::int32_t region_start;
        std::int32_t region_end;
        std::string region_chrom;
        RegionParser(region, region_chrom, region_start, region_end);

        if (region_chrom == chr){
          std::int32_t overlap = std::min(end, region_end) - std::max(start, region_start);
          overlap = std::max(0, overlap);
          if (overlap > 0){
            Fragment peak_location(chr, std::max(start, region_start), overlap, count);
            peaks.push_back(peak_location);
          }
        }
      }
    }
  }
  // sort peaks
  std::sort(peaks.begin(), peaks.end(), compare_location);
  // if peak scores have been loaded from the bed file,
  // then normalize peak scores and rescale them to 0-1
  if((count_colidx > 0) && (!noscale)) Rescale(peaks, true);
  return 0;
}

bool PeakReader::BedPeakReader(std::vector<Fragment>& peaks,
        const std::int32_t count_colidx, const std::string region, const bool noscale){
  std::ifstream infile(peakfile.c_str());
  if(infile.fail()) PrintMessageDieOnError("Input path \"" + peakfile + "\" does not exist", M_ERROR);

  std::string line;
  while (std::getline(infile, line)){
      std::string chr;
      std::int32_t start;
      std::int32_t end;
      std::int32_t length;
      float count = -1;
      
      std::stringstream linestream(line);
      std::string element;
      uint32_t elem_idx=0;
      while(!linestream.eof()){
        linestream >> element;
        if (elem_idx == 0){
          chr = element;
        }else if(elem_idx == 1){
          start = std::stol(element);
        }else if(elem_idx == 2){
          end = std::stol(element);
        }else if((count_colidx > 0) && (elem_idx == (count_colidx-1))){
          count = std::stof(element);
        }
        elem_idx++;
      }
      
      if (region.empty()){
        length = end-start;
        Fragment peak_location(chr, start, length, count);
        peaks.push_back(peak_location);
      }else{
        std::int32_t region_start;
        std::int32_t region_end;
        std::string region_chrom;
        RegionParser(region, region_chrom, region_start, region_end);

        if (region_chrom == chr){
          std::int32_t overlap = std::min(end, region_end) - std::max(start, region_start);
          overlap = std::max(0, overlap);
          if (overlap > 0){
            Fragment peak_location(chr, std::max(start, region_start), overlap, count);
            peaks.push_back(peak_location);
          }
        }
      }
  }
  // sort peaks
  std::sort(peaks.begin(), peaks.end(), compare_location);
  // if peak scores have been loaded from the bed file,
  // then normalize peak scores and rescale them to 0-1
  if((count_colidx > 0) && (!noscale)) Rescale(peaks, true);
  return 0;
}

bool PeakReader::UpdateTagCount(std::vector<Fragment>& peaks, const std::string bamfile,
            std::uint32_t* ptr_total_genome_length, float* ptr_total_tagcount, float* ptr_tagcount_in_peaks,
            const std::string region, const float frag_length){
  BamCramReader bamreader(bamfile);
  const BamHeader* bamheader = bamreader.bam_header();
  std::vector<std::string> seq_names = bamheader->seq_names();
  std::vector<uint32_t> seq_lengths = bamheader->seq_lengths();

  uint32_t total_genome_length = 0;
  std::vector<Fragment> fragments;

  if (region.empty()){
    for (int seq_index=0; seq_index<seq_names.size(); seq_index++){
      // only chr1-chr22, chrX, chrY are kept. The others are abandoned.
      if ((seq_names[seq_index].find("_") == std::string::npos) && (seq_names[seq_index] != "chrM")){
        bamreader.SetRegion(seq_names[seq_index], 0, seq_lengths[seq_index]);
        BamAlignment aln;
        float count = 0;
        while (bamreader.GetNextAlignment(aln)){
          if (aln.IsDuplicate()) {continue;} // skip the duplicated ones
          if ((!aln.IsMapped()) || aln.IsFailedQC() || aln.IsSecondary() || aln.IsSupplementary()){continue;}
          float aln_start = aln.Position();
          float aln_end = aln.GetEndPosition();
          if(aln.IsReverseStrand()){
            Fragment frag_location(seq_names[seq_index], aln_end-frag_length, frag_length, count);
            fragments.push_back(frag_location);
          }else{
            Fragment frag_location(seq_names[seq_index], aln_start, frag_length, count);
            fragments.push_back(frag_location);
          }
          count += 1;
        }

        total_genome_length += seq_lengths[seq_index];
      }
    }
  }else{
    std::int32_t region_start;
    std::int32_t region_end;
    std::string region_chrom;
    RegionParser(region, region_chrom, region_start, region_end);

    for (int seq_index=0; seq_index<seq_names.size(); seq_index++){
      if (seq_names[seq_index] == region_chrom){
        bamreader.SetRegion(seq_names[seq_index], region_start, region_end);
        BamAlignment aln;
        float count = 0;
        while (bamreader.GetNextAlignment(aln)){
          if (aln.IsDuplicate()) {continue;} // skip the duplicated ones
          if ((!aln.IsMapped()) || aln.IsFailedQC() || aln.IsSecondary() || aln.IsSupplementary()){continue;}
          float aln_start = aln.Position();
          float aln_end = aln.GetEndPosition();
          if(aln.IsReverseStrand()){
            Fragment frag_location(seq_names[seq_index], aln_end-frag_length, frag_length, count);
            fragments.push_back(frag_location);
          }else{
            Fragment frag_location(seq_names[seq_index], aln_start, frag_length, count);
            fragments.push_back(frag_location);
          }
          count += 1;
        }

        total_genome_length += (region_end-region_start);
      }
    }
  }
  *ptr_total_genome_length = total_genome_length;

  // Count number of reads that fall in peak regions
  // Naive version: (commented out)
  /*
  for (int frag_index=0; frag_index<fragments.size(); frag_index++){
    for(int peak_index = 0; peak_index < peaks.size(); peak_index++){
      if (peaks[peak_index].chromID == frags[read_index].chromID){
        if ((peaks[peak_index].start < frags[read_index].end) &&\
                (peaks[peak_index].end > frags[read_index].start)){
          n_frags_in_peak += 1;
        }
      }
    }
  }
  */
  // Faster version: (results are the same as the naive version)
  std::sort(peaks.begin(), peaks.end(), compare_location);
  std::sort(fragments.begin(), fragments.end(), compare_location);

  // reset read scores to zero
  for (int frag_index=0; frag_index<fragments.size(); frag_index++) fragments[frag_index].score = 0;
  for (int peak_index=0; peak_index<peaks.size(); peak_index++) peaks[peak_index].score = 0;

  uint32_t frag_index_start = 0;
  std::map<int, int> frag2peak;
  for(int peak_index=0; peak_index<peaks.size(); peak_index++){
    for (int frag_index=frag_index_start; frag_index<fragments.size(); frag_index++){
      if (peaks[peak_index].chrom == fragments[frag_index].chrom){
        uint32_t peak_start = peaks[peak_index].start;
        uint32_t peak_end = peaks[peak_index].start + peaks[peak_index].length;
        uint32_t frag_start = fragments[frag_index].start;
        uint32_t frag_end = fragments[frag_index].start + fragments[frag_index].length;
        if ((peak_start < frag_end) && (peak_end > frag_start)){
          float overlap = (float) (std::min(peak_end,frag_end) - std::max(peak_start, frag_start)) / (float)(frag_end-frag_start);
          fragments[frag_index].score = overlap;
          peaks[peak_index].score += overlap;
          frag2peak[frag_index] = peak_index;
        }else if(peak_start >= frag_end){
          frag_index_start = frag_index + 1;
        }else if(peak_end <= frag_start){
          break;
        }
      }else if(peaks[peak_index].chrom > fragments[frag_index].chrom){
        frag_index_start = frag_index + 1;
      }else if(peaks[peak_index].chrom < fragments[frag_index].chrom){
        break;
      }
    }
  }

  // normalize peak scores
  for(int peak_index=0; peak_index<peaks.size(); peak_index++){
    peaks[peak_index].score /= ((float) peaks[peak_index].length);
  }
  Rescale(peaks, true);

  // total num of fragments in peaks
  float n_frags_in_peak = 0;
  for (int frag_index=0; frag_index<fragments.size(); frag_index++){
    int peak_index = frag2peak[frag_index];
    n_frags_in_peak += (fragments[frag_index].score * peaks[peak_index].score);
  }

  *ptr_tagcount_in_peaks = n_frags_in_peak;
  *ptr_total_tagcount = (float) fragments.size();
  return 0;
}

void PeakReader::Rescale(std::vector<Fragment>& peaks, bool rm_outliers) {
  // Find max, median
  std::vector<float> scores;
  float max_score = 0;
  for(int peak_index=0; peak_index<peaks.size(); peak_index++){
    scores.push_back(peaks[peak_index].score);
    if (peaks[peak_index].score > max_score) max_score = peaks[peak_index].score;
  }
  float threshold = max_score;
  if (rm_outliers) {
    size_t n = scores.size()/2;
    nth_element(scores.begin(), scores.begin()+n, scores.end());
    float median_score = scores[n];
    float max_possible = median_score*2;
    if (max_possible<max_score) {
      threshold = max_possible;
    }
  }
  for(int peak_index=0; peak_index<peaks.size(); peak_index++) {
    if (peaks[peak_index].score >= threshold) {
      peaks[peak_index].score = 1;
    } else {
      peaks[peak_index].score /= threshold;
    }
  }
}

bool PeakReader::compare_location(Fragment a, Fragment b){
  if (a.chrom != b.chrom){
    return (a.chrom < b.chrom);
  }else if (a.start != b.start){
    return (a.start < b.start);
  }else if (a.length != b.length){
    return (a.length < b.length);
  }else{
    return (a.score < b.score);
  }
}

void PeakReader::RegionParser(const std::string region,
        std::string& chromID, std::int32_t& start, std::int32_t& end){
  std::stringstream ss(region);

  // parse chromID and the rest
  std::string split;
  std::vector<std::string> parts;
  while(std::getline(ss, split, ':')) parts.push_back(split);
  if (parts.size() != 2) PrintMessageDieOnError("Improper region input format should be chrom:start-end", M_ERROR);
  chromID = parts[0];
  std::string start_end = parts[1];

  // reset
  parts.clear();
  ss.str(start_end);
  ss.clear();

  // parse start and end
  while(getline(ss, split, '-')) parts.push_back(split);
  if (parts.size() != 2) PrintMessageDieOnError("Improper region input format should be chrom:start-end", M_ERROR);
  start = stoi(parts[0]);
  end = stoi(parts[1]);
}

/* END OF FILE */
