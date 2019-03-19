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
  chroms->clear();
  int nseqs = faidx_nseq(refindex);
  for (int i=0; i<nseqs; i++) {
    chroms->push_back(string(faidx_iseq(refindex, i)));
  }
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
  vector<string> chroms;
  if (!GetChroms(&chroms)) {
    return false;
  }
  for (size_t i=0; i<chroms.size(); i++) {
    int length = faidx_seq_len(refindex, chroms[i].c_str());
    (*chromLengths)[chroms[i]] = length;
  }
  return true;
}

RefGenome::~RefGenome() {
  fai_destroy(refindex);
}

int32_t RefGenome::GetGenomeLength() {
  std::map<std::string, int> chromlens;
  if (!GetLengths(&chromlens)) {
    PrintMessageDieOnError("Error getting chrom lengths", M_ERROR);
  }
  int32_t gsize = 0;
  for (std::map<std::string, int>::iterator it = chromlens.begin(); it != chromlens.end(); it++) {
    gsize += it->second;
  }
  return gsize;
}
