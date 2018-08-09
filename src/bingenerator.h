#ifndef SRC_BINGENERATOR_H__
#define SRC_BINGENERATOR_H__

#include "src/options.h"

using namespace std;

class GenomeBin {
 public:
  std::string chrom;
  int32_t start, end;

  GenomeBin(std::string chrom_, int32_t start_, int32_t end_) {
    chrom = chrom_;
    start = start_;
    end = end_;
  }

  ~GenomeBin() {}
};

class BinGenerator {
  /*
    This class manages binning the genome
    If options.region is set, only get bins from that region
    Otherwise, bin the entire genome, getting the size from options.reffa
    Generate bins of size options.binsize

    Before implementing, take a look at bedtools makewindows
    Also, see ref_genome.h for reference genome class you will probably have to use here
   */
 public:
  BinGenerator(const Options& options);
  virtual ~BinGenerator();

  /* Set the next bin */
  bool GotoNextBin();

  /* Return the current bin */
  const GenomeBin GetCurrentBin();

 private:
  GenomeBin* currentBin;
  int32_t reg_end;
  int binsize;
};

#endif  // SRC_BINGENERATOR_H__
