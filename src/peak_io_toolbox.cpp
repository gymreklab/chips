// Peak Loader Tool Box
//

#include "src/peak_io_toolbox.h"

PeakReader::PeakReader(const std::string& _peakfile){
  peakfile = _peakfile;
}

bool PeakReader::HomerPeakReader(std::vector<Fragment>& peaks, const std::int32_t count_colidx){
  std::ifstream infile(peakfile.c_str());
  std::string line;
  while (std::getline(infile, line)){
    if (line.find('#') != 0){
      std::string chr;
      std::uint32_t start;
      std::uint32_t end;
      std::uint32_t length;
      float count = -1;
      
      std::stringstream linestream(line);
      std::string element;
      uint32_t elem_idx=0;
      while(!linestream.eof()){
        linestream >> element;
        if (elem_idx == 1){
          chr = element;
        }else if(elem_idx == 2){
          start = std::stoul(element);
        }else if(elem_idx == 3){
          end = std::stoul(element);
        }else if((count_colidx != -1) && (elem_idx == count_colidx)){
          count = std::stof(element);
        }
        elem_idx++;
      }
      
      length = end-start;
      Fragment peak_location(chr, start, length, count);
      peaks.push_back(peak_location);
    }
  }
  return 0;
}

bool PeakReader::TestPeakReader(std::vector<Fragment>& peaks, const std::int32_t count_colidx){
  std::ifstream infile(peakfile.c_str());
  std::string line;
  while (std::getline(infile, line)){
      std::string chr;
      std::uint32_t start;
      std::uint32_t end;
      std::uint32_t length;
      float count = -1;
      
      std::stringstream linestream(line);
      std::string element;
      uint32_t elem_idx=0;
      while(!linestream.eof()){
        linestream >> element;
        if (elem_idx == 0){
          chr = element;
        }else if(elem_idx == 1){
          start = std::stoul(element);
        }else if(elem_idx == 2){
          end = std::stoul(element);
        }else if((count_colidx != -1) && (elem_idx == count_colidx)){
          count = std::stof(element);
          //std::cout<< count<<std::endl;
        }
        elem_idx++;
      }
      
      length = end-start;
      Fragment peak_location(chr, start, length, count);
      peaks.push_back(peak_location);
  }
  return 0;
}

bool PeakReader::UpdateTagCount(std::vector<Fragment>& peaks, const std::string bamfile,
            uint32_t* ptr_total_genome_length, uint32_t* ptr_total_tagcount, uint32_t* ptr_tagcount_in_peaks){
  BamCramReader bamreader(bamfile);
  const BamHeader* bamheader = bamreader.bam_header();
  std::vector<std::string> seq_names = bamheader->seq_names();
  std::vector<uint32_t> seq_lengths = bamheader->seq_lengths();

  uint32_t total_genome_length = 0;
  std::vector<Fragment> reads;
  for (int seq_index=0; seq_index<seq_names.size(); seq_index++){
    // only chr1-chr22, chrX, chrY are kept. The others are abandoned.
    if ((seq_names[seq_index].find("_") == std::string::npos) && (seq_names[seq_index] != "chrM")){
      bamreader.SetRegion(seq_names[seq_index], 0, seq_lengths[seq_index]);
      BamAlignment aln;
      while (bamreader.GetNextAlignment(aln)){
        Fragment read_location(seq_names[seq_index], aln.Position(), aln.Length(), 0);
        reads.push_back(read_location);
      }

      total_genome_length += seq_lengths[seq_index];
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
/* END OF FILE */
