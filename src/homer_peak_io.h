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

/* data structure */
struct seq_loc{
    std::string chromID;
    uint32_t start;
    uint32_t end;
};

/* Peak File Reader */
class HomerPeakReader{
  public:
    HomerPeakReader(const std::string& peakfile, std::vector<seq_loc>& peaks);
};

#endif
