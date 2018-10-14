#include "src/sequencer.h"
#include <iostream>
#include <sstream>

const std::map<char, char> Sequencer::NucleotideMap = {
                    {'A', 'T'}, {'C', 'G'}, {'T', 'A'}, {'G', 'C'}, {'N', 'N'},
                    {'a', 't'}, {'c', 'g'}, {'t', 'a'}, {'g', 'c'}, {'n', 'n'}};

Sequencer::Sequencer(const Options& options) {
  ref_genome = new RefGenome(options.reffa);
  paired = options.paired;
  outprefix = options.outprefix;
  readlen = options.readlen;
}

void Sequencer::Sequence(const std::vector<Fragment>& input_fragments, int& fastq_index, int thread_index, int copy_index) {
  std::string frag_seq;
  std::string read_seq;
  std::string read_seq_rc;
  std::vector<std::string> read_pair;
  std::vector<std::string> reads_1;
  std::vector<std::string> reads_2;
  std::vector<std::string> chroms;
  std::vector<std::string> ids;
  std::vector<int> starts_1;
  std::vector<int> starts_2;

  for (int frag_index=0; frag_index<input_fragments.size(); frag_index++){
    if(ref_genome->GetSequence(input_fragments[frag_index].chrom,
            input_fragments[frag_index].start,
            input_fragments[frag_index].start+input_fragments[frag_index].length,
            &frag_seq)){

      // generate reads from both strands
      read_pair.clear();
      if (Fragment2Read(frag_seq, read_seq) &&
            Fragment2Read(ReverseComplement(frag_seq), read_seq_rc)){
        read_pair.push_back(read_seq);
        read_pair.push_back(read_seq_rc);
	std::stringstream ss;
	ss << input_fragments[frag_index].chrom << ":"
	   << input_fragments[frag_index].start << ":"
	   << input_fragments[frag_index].length;
	ids.push_back(ss.str());
	std::random_shuffle(read_pair.begin(), read_pair.end());
      }else{
        continue;
      }
      chroms.push_back(input_fragments[frag_index].chrom);
      starts_1.push_back(input_fragments[frag_index].start);
      starts_2.push_back(input_fragments[frag_index].start+input_fragments[frag_index].length-readlen+1);
      reads_1.push_back(read_pair[0]);
      reads_2.push_back(read_pair[1]);
    }
  }

  // save into file
  if (paired){
    int temp = fastq_index;
    save_into_fastq(reads_1, ids, outprefix + "_"+std::to_string(thread_index)+"_1.fastq", temp, copy_index);
    save_into_fastq(reads_2, ids, outprefix + "_"+std::to_string(thread_index)+"_2.fastq", fastq_index, copy_index);
    //    save_into_sam(reads_1, reads_2, chroms, starts_1, starts_2, outprefix + "/reads_"+std::to_string(thread_index)+".sam");
  }else{
    save_into_fastq(reads_1, ids, outprefix + "_"+std::to_string(thread_index)+".fastq", fastq_index, copy_index);
  }
}

bool Sequencer::Fragment2Read(const std::string frag, std::string& read){
  try{
    read = frag.substr(0, readlen);
    if (frag.length() < readlen){ return false;}
    return true;
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
