#include "src/sequencer.h"

const std::map<char, char> Sequencer::NucleotideMap = {
                    {'A', 'T'}, {'C', 'G'}, {'T', 'A'}, {'G', 'C'}, {'N', 'N'},
                    {'a', 't'}, {'c', 'g'}, {'t', 'a'}, {'g', 'c'}, {'n', 'n'}};

Sequencer::Sequencer(const Options& options) {
  ref_genome = new RefGenome(options.reffa);
  paired = options.paired;
  outprefix = options.outprefix;
  readlen = options.readlen;
}

void Sequencer::Sequence(const std::vector<Fragment>& input_fragments, int* fastq_index) {
  std::string frag_seq;
  std::string read_seq;
  std::string read_seq_rc;
  std::vector<std::string> read_pair;
  std::vector<std::string> reads_1;
  std::vector<std::string> reads_2;

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
        std::random_shuffle(read_pair.begin(), read_pair.end());
      }else{
        continue;
      }

      reads_1.push_back(read_pair[0]);
      reads_2.push_back(read_pair[1]);
    }
  }

  // save into file
  if (paired){
    int temp = *fastq_index;
    save_into_fastq(reads_1, outprefix + "/reads_1.fastq", &temp);
    save_into_fastq(reads_2, outprefix + "/reads_2.fastq", fastq_index);
  }else{
    save_into_fastq(reads_1, outprefix + "/reads.fastq", fastq_index);
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

// TODO update this function so the read names are unique
bool Sequencer::save_into_fastq(const std::vector<std::string> reads, const std::string ofilename, int* fastq_index){
  std::ofstream ofile(ofilename, std::ofstream::out | std::ofstream::app);
  for (int read_index=0; read_index<reads.size(); read_index++){
    ofile << "@SIM:" << read_index + *fastq_index <<"\n";
    ofile << reads[read_index] << "\n";
    ofile << "+\n";
    ofile << std::string(readlen, '~') << "\n";
  }
  ofile.close();
  *fastq_index += reads.size();
  return true;
}

Sequencer::~Sequencer() {
  delete ref_genome;
}
