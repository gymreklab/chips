#ifndef SRC_SEQUENCER_H__
#define SRC_SEQUENCER_H__

#include "src/fragment.h"
#include "src/options.h"
#include "src/ref_genome.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>

class Sequencer {
 public:
  Sequencer(const Options& options);
  virtual ~Sequencer();

  void Sequence(const std::vector<Fragment>& input_fragments, \
                    int& fastq_index, int thread_index, int copy_index);
 private:
  RefGenome* ref_genome;
  bool paired;
  std::string outprefix;
  int readlen;
  static const std::map<char, char> NucleotideMap;

  bool Fragment2Read(const std::string frag, std::string& read);
  std::string ReverseComplement(const std::string seq);
  bool save_into_fastq(const std::vector<std::string> reads,
		       const std::vector<std::string> ids,
		       const std::string ofilename,
                        int& fastq_index, int copy_index);
  /*
  bool save_into_sam(const std::vector<std::string> reads1,
		     const std::vector<std::string> reads2,
		     const std::vector<std::string> chroms,
		     const std::vector<int> starts1,
		     const std::vector<int> starts2,
		     const std::string ofilename);
  */
};

#endif  // SRC_SEQUENCER_H__
