#include "sequencer.h"
#include "common.h"
#include <iostream>
#include <sstream>
#include <random>
#include <algorithm>
#include <iterator>

const std::map<char, char> Sequencer::NucleotideMap = {
                    {'A', 'T'}, {'C', 'G'}, {'T', 'A'}, {'G', 'C'}, {'N', 'N'},
                    {'a', 't'}, {'c', 'g'}, {'t', 'a'}, {'g', 'c'}, {'n', 'n'}};

const char Sequencer::NucleotideTypesUpper[] = {'A','T','C','G'};
const char Sequencer::NucleotideTypesLower[] = {'a','t','c','g'};

const std::map<char, std::vector<char> > Sequencer::SubMap{
                    {'A', {'T','C','G'}},
                    {'T', {'A','C','G'}},
                    {'C', {'A','T','G'}},
                    {'G', {'A','T','C'}},
                    {'a', {'t','c','g'}},
                    {'t', {'a','c','g'}},
                    {'c', {'a','t','g'}},
                    {'g', {'a','t','c'}}};

Sequencer::Sequencer(const Options& options) {
  ref_genome = new RefGenome(options.reffa);
  paired = options.paired;
  outprefix = options.outprefix;
  readlen = options.readlen;
  pcr_rate = options.pcr_rate;

  sequencer_type = options.sequencer_type;
  if (sequencer_type == ""){
    sub_rate = options.sub_rate;
    del_rate = options.del_rate;
    ins_rate = options.ins_rate;
  }else if(sequencer_type == "HiSeq"){
    //https://academic.oup.com/bib/article/17/1/154/1742474
    sub_rate = 2.640*1e-3;
    del_rate = 2.433*1e-4;
    ins_rate = 1.833*1e-4;
  }
}


void Sequencer::Sequence(const std::vector<Fragment>& input_fragments,
			 const int& numreads,
			 int& fastq_index, int thread_index, int copy_index, std::mt19937& rng) {
  std::string frag_seq;
  std::string read_seq;
  std::string read_seq_rc;
  std::vector<std::string> read_pair;
  std::vector<std::string> reads_1;
  std::vector<std::string> reads_2;
  std::vector<std::string> ids;

  // Sample from fragments w/o replacement (by shuffling first)
  // If needed, go through the fragments multiple times
  int total_reads_sequenced = 0;
  std::vector<size_t> frag_indices;
  size_t frag_index;
  for (size_t frag_index=0; frag_index<input_fragments.size(); frag_index++) {
    frag_indices.push_back(frag_index);
  }

  std::stringstream ss;
  ss << "Sequencing " << numreads << " total reads from " << input_fragments.size() << " fragments for copy " << copy_index;
  int dupnum = 0;
  //  PrintMessageDieOnError(ss.str(), M_PROGRESS);
  while (true) {
    std::shuffle(frag_indices.begin(), frag_indices.end(), rng);
    for (size_t fg=0; fg<frag_indices.size(); fg++) {
      frag_index = frag_indices[fg];
      if(!ref_genome->GetSequence(input_fragments[frag_index].chrom,
				 input_fragments[frag_index].start,
				 input_fragments[frag_index].start+input_fragments[frag_index].length,
				  &frag_seq)) {
	continue;
      }

      dupnum = 0;
      // generate reads from both strands
      read_pair.clear();
      if (Fragment2Read(frag_seq, read_seq, rng) &&
	  Fragment2Read(ReverseComplement(frag_seq), read_seq_rc, rng)){
	read_pair.push_back(read_seq);
	read_pair.push_back(read_seq_rc);
	std::stringstream ss;
	ss << input_fragments[frag_index].chrom << ":"
	   << input_fragments[frag_index].start << ":"
	   << input_fragments[frag_index].length << ":";
	std::shuffle(read_pair.begin(), read_pair.end(), rng);
	// Add the first read
	reads_1.push_back(read_pair[0]);
	reads_2.push_back(read_pair[1]);
	ids.push_back(ss.str()+std::to_string(dupnum));
	total_reads_sequenced += 1;
	dupnum += 1;
      } else {
	continue;
      }

      // PCR 
      while (true) {
        if (total_reads_sequenced >= numreads) break;
        if ( ((float) rng()/(float) rng.max()) < pcr_rate) break;
        ids.push_back(ids[ids.size()-1]+std::to_string(dupnum));
        reads_1.push_back(read_pair[0]);
        reads_2.push_back(read_pair[1]);
        total_reads_sequenced += 1;
	dupnum += 1;
      } 

      if (total_reads_sequenced >= numreads) {
	    break;
      }
    }
    if (total_reads_sequenced >= numreads) {
      break;
    }
  }
    
  // save into file
  if (paired) {
    int temp = fastq_index;
    save_into_fastq(reads_1, ids, outprefix + "_"+std::to_string(thread_index)+"_1.fastq", temp, copy_index);
    save_into_fastq(reads_2, ids, outprefix + "_"+std::to_string(thread_index)+"_2.fastq", fastq_index, copy_index);
    //    save_into_sam(reads_1, reads_2, chroms, starts_1, starts_2, outprefix + "/reads_"+std::to_string(thread_index)+".sam");
  } else{
    save_into_fastq(reads_1, ids, outprefix + "_"+std::to_string(thread_index)+".fastq", fastq_index, copy_index);
  }
}

bool Sequencer::Fragment2Read(const std::string frag, std::string& read, std::mt19937& rng){
  try{
    //read = frag.substr(0, readlen);
    //if (frag.length() < readlen){ return false;}
    read = "";
    double dice;
    int elem_index = 0;
    while (read.size() < readlen){
      dice = (((float) rng()/(float) rng.max()));
      if (dice <= ins_rate){
        // randomly insert a nucleotide
        int ins_index = rng() % 4;
        read += NucleotideTypesUpper[ins_index];
      }else if (dice <= (ins_rate + del_rate)){
        // skip this nucleotide
        if (elem_index < frag.size()){
          elem_index += 1;
          continue;
        }else{
          break;
        }
      }else if (dice <= (ins_rate + del_rate + sub_rate)){
        // substitude this nucleotide with another
        if (elem_index < frag.size()){
          char nuc_to_mut = frag[elem_index];
          if ((nuc_to_mut == 'N') || (nuc_to_mut == 'n')){
            read += frag[elem_index];
          }else{
            int sub_index = rng() % 3;
            read += (SubMap.at(nuc_to_mut))[sub_index];
          }
          elem_index += 1;
        }else{
          break;
        }
      }else{
        if (elem_index < frag.size()){
          // correctly sequenced
          read += frag[elem_index];
          elem_index += 1;
        }else{
          break;
        }
      }
    }

    if (read.size() < readlen){
      // fill up the reads with "N"s if the fragment length is
      // shorter than the read length
      int N_total = readlen-read.size();
      for (int add_n=0; add_n<N_total; add_n++) read += 'N';
      return true;
    }else if (read.size() == readlen){
      return true;
    }else {
      std::cerr << "Unexpected read size: " << read.size() << std::endl;
      return false;
    }
  } catch (const char* msg){
    std::cerr << msg << " in Sequencer::Fragment2Read!"<< std::endl;
    return false;
  }
}

std::string Sequencer::ReverseComplement(const std::string seq){
  std::string rc_seq;
  for (int i=(seq.length()-1); i>=0; i--){
    rc_seq += NucleotideMap.at(seq[i]);
  }
  return rc_seq; 
}

/*
bool Sequencer::save_into_sam(const std::vector<std::string> reads1,
			      const std::vector<std::string> reads2,
			      const std::vector<std::string> chroms,
			      const std::vector<int> starts1,
			      const std::vector<int> starts2,
			      const std::string ofilename) {
  std::ofstream ofile(ofilename, std::ofstream::out | std::ofstream::app);
  // ofile << "@HD\tVN:1.0\tSO:coordinate\n";
  //ofile << "@SQ\tSN:chr19\tLN:59128983\n";
  for (int read_index=0; read_index<reads1.size(); read_index++){
    ofile << "asimon:" << chroms[read_index] << ":" << starts1[read_index] << "\t"
	  << 67 << "\t"
	  << chroms[read_index] << "\t"
	  << starts1[read_index] << "\t"
	  << 255 << "\t"
	  << readlen << "M" << "\t"
	  << "=" << "\t"
	  << starts2[read_index] << "\t"
	  << starts2[read_index]+readlen-starts1[read_index] << "\t"
	  << reads1[read_index] << "\t"
	  << std::string(readlen, '@') << "\n";
    ofile << "asimon:" << chroms[read_index] << ":" << starts1[read_index] << "\t"
	  << 147 << "\t"
	  << chroms[read_index] << "\t"
	  << starts2[read_index] << "\t"
	  << 255 << "\t"
	  << readlen << "M" << "\t"
	  << "=" << "\t"
	  << starts1[read_index] << "\t"
	  << starts2[read_index]+readlen-starts1[read_index] << "\t"
	  << ReverseComplement(reads2[read_index]) << "\t"
	  << std::string(readlen, '@') << "\n";
  }
  ofile.close();
  return true;
  }*/

bool Sequencer::save_into_fastq(const std::vector<std::string> reads, 
				const std::vector<std::string> ids,
				const std::string ofilename,
				int& fastq_index, int copy_index){
  std::ofstream ofile(ofilename, std::ofstream::out | std::ofstream::app);
  for (int read_index=0; read_index<reads.size(); read_index++){
    ofile << "@SIM:" << ids[read_index] << ":" << copy_index<<":"<< read_index + fastq_index <<"\n";
    ofile << reads[read_index] << "\n";
    ofile << "+\n";
    ofile << std::string(readlen, '~') << "\n";
  }
  ofile.close();
  fastq_index += reads.size();
  return true;
}

Sequencer::~Sequencer() {
  delete ref_genome;
}
