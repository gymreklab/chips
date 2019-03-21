/*
 * Created by An Zheng
 * Modified by An Zheng
 * Last modified in July 29th, 2018
 * */

#ifndef PEAK_LOAD_H_
#define PEAK_LOAD_H_

#include <iostream>
#include <inttypes.h>
#include <stdbool.h>
#include <assert.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <map>

#include "src/peak_io_toolbox.h"
#include "src/fragment.h"

class PeakLoader{
  public:
    PeakLoader(const std::string _peakfile, const std::string _peakfileType="",
                        const std::string _bamfile="", const std::int32_t _count_colidx=-1);
    bool Load(std::vector<Fragment>& peaks, const std::string region="", const float frag_length=0, const bool noscale=false);
    std::uint32_t total_genome_length;
    float total_tagcount;
    float tagcount_in_peaks;
  private:
    static const std::map<std::string, int> peakfileTypeList;
    std::string peakfile;
    std::string peakfileType;
    std::string bamfile;
    std::uint32_t count_colidx;
};
#endif
