/*
 * Created by An Zheng
 * Modified by An Zheng
 * Last modified in July 29th, 2018
 * */

#ifndef PEAK_IO_H_
#define PEAK_IO_H_

#include <iostream>
#include <inttypes.h>
#include <stdbool.h>
#include <assert.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include "src/fragment.h"
#include "src/bam_io.h"

/* Peak File Reader */
class PeakReader{
  public:
    PeakReader(const std::string& peakfile);
    bool HomerPeakReader(std::vector<Fragment>& peaks, const std::int32_t count_colidx);
    bool TestPeakReader(std::vector<Fragment>& peaks, const std::int32_t count_colidx);
    bool UpdateTagCount(std::vector<Fragment>& peaks, const std::string bamfile,
            uint32_t* ptr_total_genome_length, uint32_t* ptr_total_tagcount, uint32_t* ptr_tagcount_in_peaks);
  private:
    std::string peakfile;
    static bool compare_location(Fragment a, Fragment b);
};


#endif
