#include <stdlib.h>

#include <iostream>
#include <sstream>
#include <algorithm>

#include "src/common.h"
#include "src/ref_genome.h"

using namespace std;

RefGenome::RefGenome(const std::string& _reffa) {
  // Check if file exists
  if (!file_exists(_reffa)) {
    PrintMessageDieOnError("FASTA file " + _reffa + " does not exist", M_ERROR);
  }

  // Check for index
  if (!file_exists(_reffa + ".fai")) {
    PrintMessageDieOnError("No index for FASTA file " + _reffa, M_ERROR);
  }

  // Load index
  refindex = fai_load(_reffa.c_str());
  if (refindex == NULL) {
    PrintMessageDieOnError("Failed to load FASTA index", M_ERROR);
  }
}

bool RefGenome::GetSequence(const std::string& _chrom,
			    const int32_t& _start,
			    const int32_t& _end,
			    std::string* seq) {
  int length;
  char* result = faidx_fetch_seq(refindex, _chrom.c_str(), _start, _end, &length);
  if (result == NULL) {
    stringstream ss;
    ss << "Error fetching reference sequence for " << _chrom << ":" << _start;
    PrintMessageDieOnError(ss.str(), M_ERROR);
  }
  seq->assign(result, length);
  std::transform(seq->begin(), seq->end(), seq->begin(), ::tolower);
  free((void *)result);
  return true;
}

/*
   Inputs:
    - chroms: vector pointer where all of the chromosomes will be stored.
 
   Outputs: 
    - bool: return true if getting the chromosomes was successful.
            False otherwise.
  
   Get a list of chromosomes from the reference genome.
*/

bool RefGenome::GetChroms(vector<string>* chroms) {
  // TODO Update this function.
  // dummy chroms list for testing bingenerator
  *chroms = {"chr_1", "chr_2", "chr_3", "chr_4", "chr_5", "chr_6", "chr_7",
             "chr_8", "chr_9", "chr_10", "chr_11", "chr_12", "chr_13",
             "chr_14", "chr_15"};
  return true;
}

/*
   Inputs:
    - chromLengths: map pointer that stores key value pairs of
                    {Chromosome, total length}.
 
   Outputs:
    - bool: return true if getting the chromosome lengths was successful.
            False otherwise.
  
   Get a map of {chromosomes, total lengths} from the reference genome.
*/
bool RefGenome::GetLengths(map<string, int>* chromLengths) {
  // TODO Update this function.
  // dummy chroms, lengths map for testing bingenerator
  *chromLengths = {{"chr_1", 248956422},{"chr_2", 242193529},{"chr_3", 198295559},
                  {"chr_4", 190214555},{"chr_5", 181538259},{"chr_6", 170805979},
                  {"chr_7", 159345973},{"chr_8", 145138636},{"chr_9", 138394717},
                  {"chr_10", 133797422},{"chr_11", 135086622},{"chr_12", 133275309},
                  {"chr_13", 114364328},{"chr_14", 107043718},{"chr_15", 101991189}};
  return true;
}

RefGenome::~RefGenome() {
  fai_destroy(refindex);
}
