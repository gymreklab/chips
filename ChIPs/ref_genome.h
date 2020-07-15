#ifndef SRC_REF_GENOME_H__
#define SRC_REF_GENOME_H__

#include "htslib/faidx.h"

#include <stdint.h>
#include <unistd.h>

#include <map>
#include <vector>
#include <string>

class RefGenome {
 public:
  RefGenome(const std::string& _reffa);
  virtual ~RefGenome();

  bool GetSequence(const std::string& _chrom,
		   const int32_t& _start,
		   const int32_t& _end,
		   std::string* seq);

  bool GetChroms(std::vector<std::string>* chroms);

  bool GetLengths(std::map<std::string, int>* chromLengths);

  int64_t GetGenomeLength();
 private:
  bool file_exists(std::string path) const {
    return (access(path.c_str(), F_OK) != -1);
  }

  faidx_t* refindex;
};

#endif  // SRC_REF_GENOME_H__
