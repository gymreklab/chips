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
#include "fragment.h"
#include "bam_io.h"

/* Peak File Reader */
class PeakReader{
  public:
    PeakReader(const std::string& peakfile);
    bool HomerPeakReader(std::vector<Fragment>& peaks,
			 const std::int32_t count_colidx, const std::string region, const bool noscale, const bool scale_outliers);
    bool BedPeakReader(std::vector<Fragment>& peaks, const std::int32_t count_colidx, const std::string region, const bool noscale, const bool scale_outliers);
    bool EmptyPeakReader();
    bool UpdateTagCount(std::vector<Fragment>& peaks, const std::string bamfile,
			std::uint32_t* ptr_total_genome_length, float* ptr_total_tagcount,
			float* ptr_tagcount_in_peaks,const std::string region, const float frag_length,
			const bool noscale, const bool scale_outliers);

  private:
    std::string peakfile;
    static bool compare_location(Fragment a, Fragment b);
    void Rescale(std::vector<Fragment>& peaks, bool rm_outliers);
};


#endif
