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
    bool HomerPeakReader(std::vector<Fragment>& peaks,
            const std::int32_t count_colidx, const std::string region);
    bool BedPeakReader(std::vector<Fragment>& peaks, const std::int32_t count_colidx, const std::string region);
    bool TestPeakReader(std::vector<Fragment>& peaks, const std::int32_t count_colidx, const std::string region);
    bool UpdateTagCount(std::vector<Fragment>& peaks, const std::string bamfile,
            std::uint32_t* ptr_total_genome_length, float* ptr_total_tagcount,
            float* ptr_tagcount_in_peaks,const std::string region, const float frag_length);
  private:
    std::string peakfile;
    static bool compare_location(Fragment a, Fragment b);
    static void RegionParser(const std::string region, std::string& chromID, std::int32_t& start, std::int32_t& end);
};


#endif
