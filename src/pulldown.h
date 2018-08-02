#ifndef SRC_PULLDOWN_H__
#define SRC_PULLDOWN_H__

#include "src/bingenerator.h"
#include "src/fragment.h"
#include "src/options.h"

#include <vector>

using namespace std;

class Pulldown {
 public:
  Pulldown(const Options& options, const GenomeBin& gbin);
  virtual ~Pulldown();

  void Perform(vector<Fragment>* output_fragments);

 private:
  std::string chrom;
  int32_t start;
  int32_t end;
  int numcopies;
  float gamma_alpha, gamma_beta;
  float ratio_beta;
};

#endif  // SRC_PULLDOWN_H__
