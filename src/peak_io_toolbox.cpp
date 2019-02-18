// Peak Loader Tool Box
//

#include "src/peak_io_toolbox.h"

PeakReader::PeakReader(const std::string& _peakfile){
  peakfile = _peakfile;
}

bool PeakReader::HomerPeakReader(std::vector<Fragment>& peaks,
        const std::int32_t count_colidx, const std::string region){
  std::ifstream infile(peakfile.c_str());
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
  std::sort(peaks.begin(), peaks.end(), compare_location);
  return 0;
}

bool PeakReader::BedPeakReader(std::vector<Fragment>& peaks,
        const std::int32_t count_colidx, const std::string region){
  std::ifstream infile(peakfile.c_str());
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
          //std::cout<< count<<std::endl;
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
  std::sort(peaks.begin(), peaks.end(), compare_location);
  return 0;
}

bool PeakReader::TestPeakReader(std::vector<Fragment>& peaks,
        const std::int32_t count_colidx, const std::string region){
  std::ifstream infile(peakfile.c_str());
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
        }else if((count_colidx > 0) && (elem_idx == count_colidx-1)){
          count = std::stof(element);
          //std::cout<< count<<std::endl;
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
  std::sort(peaks.begin(), peaks.end(), compare_location);
  return 0;
}

bool PeakReader::UpdateTagCount(std::vector<Fragment>& peaks, const std::string bamfile,
            uint32_t* ptr_total_genome_length, uint32_t* ptr_total_tagcount, uint32_t* ptr_tagcount_in_peaks,
            const std::string region){
  BamCramReader bamreader(bamfile);
  const BamHeader* bamheader = bamreader.bam_header();
  std::vector<std::string> seq_names = bamheader->seq_names();
  std::vector<uint32_t> seq_lengths = bamheader->seq_lengths();

  uint32_t total_genome_length = 0;
  std::vector<Fragment> reads;

  if (region.empty()){
    for (int seq_index=0; seq_index<seq_names.size(); seq_index++){
      // only chr1-chr22, chrX, chrY are kept. The others are abandoned.
      if ((seq_names[seq_index].find("_") == std::string::npos) && (seq_names[seq_index] != "chrM")){
        bamreader.SetRegion(seq_names[seq_index], 0, seq_lengths[seq_index]);
        BamAlignment aln;
        while (bamreader.GetNextAlignment(aln)){
          if (aln.IsDuplicate()) {continue;} // skip the duplicated ones
          if ( (!aln.IsMapped()) || aln.IsSecondary()){continue;}
          Fragment read_location(seq_names[seq_index], aln.Position(), aln.Length(), 0);
          reads.push_back(read_location);
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
        while (bamreader.GetNextAlignment(aln)){
          if (aln.IsDuplicate()) {continue;} // skip the duplicated ones
          if ( (!aln.IsMapped()) || aln.IsSecondary()){continue;}
          Fragment read_location(seq_names[seq_index], aln.Position(), aln.Length(), 0);
          reads.push_back(read_location);
        }

        total_genome_length += (region_end-region_start);
      }
    }
  }
  *ptr_total_genome_length = total_genome_length;

  // Count number of reads that fall in peak regions
  // Naive version: (commented out)
  /*
  for (int read_index=0; read_index<reads.size(); read_index++){
    for(int peak_index = 0; peak_index < peaks.size(); peak_index++){
      if (peaks[peak_index].chromID == reads[read_index].chromID){
        if ((peaks[peak_index].start < reads[read_index].end) &&\
                (peaks[peak_index].end > reads[read_index].start)){
          n_reads_in_peak += 1;
        }
      }
    }
  }
  */
  // Faster version: (results are the same as the naive version)
  std::sort(peaks.begin(), peaks.end(), compare_location);
  std::sort(reads.begin(), reads.end(), compare_location);

  uint32_t read_index_start = 0;
  for(int peak_index=0; peak_index<peaks.size(); peak_index++){
    for (int read_index=read_index_start; read_index<reads.size(); read_index++){
      if (peaks[peak_index].chrom == reads[read_index].chrom){
        uint32_t peak_start = peaks[peak_index].start;
        uint32_t peak_end = peaks[peak_index].start + peaks[peak_index].length;
        uint32_t read_start = reads[read_index].start;
        uint32_t read_end = reads[read_index].start + reads[read_index].length;
        if ((peak_start < read_end) && (peak_end > read_start)){
          reads[read_index].score = 1;
          peaks[peak_index].score += 1;
          //n_reads_in_peak += 1;
          //break;
        }else if(peak_start >= read_end){
          read_index_start = read_index + 1;
        }else if(peak_end <= read_start){
          break;
        }
      }else if(peaks[peak_index].chrom > reads[read_index].chrom){
        read_index_start = read_index + 1;
      }else if(peaks[peak_index].chrom < reads[read_index].chrom){
        break;
      }
    }
  }

  uint32_t n_reads_in_peak = 0;
  for (int read_index=0; read_index<reads.size(); read_index++){
    n_reads_in_peak += float(reads[read_index].score);
  }

  *ptr_tagcount_in_peaks = n_reads_in_peak;
  *ptr_total_tagcount = reads.size();
  return 0;
}

bool PeakReader::compare_location(Fragment a, Fragment b){
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
