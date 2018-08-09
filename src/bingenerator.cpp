#include "src/bingenerator.h"
#include "src/common.h"

#include <sstream>
#include <iostream>

/*
  Constructor for BinGenerator
  
  Set up binning over a region, or over the entire genome
  Initialize currentBin to the first bin
 */
BinGenerator::BinGenerator(const Options& options) {
  // Find the appropriate step size good way to check if in region or not
  // TODO Region format chrom:start-end will be divided into bins of size options.binsize
  binsize = options.binsize;
  std::string chrom;
  int32_t start;
  int32_t end;
  
  if (!optons.region || option.region == "genome-wide")
  {
    // default is genome wide TODO figure out if its null or "genome-wide" as default
    // TODO How to handle traversing through whole genome?
    //      - check if last chromosome and 
  }
  else
  {
    std::vector<std::string> parts;        // vector to store split strings
    std::istringstream ss(options.region); // read in region
    std::string split;                     // used to store split strings

    // get chrom and locations from string
    while(std::getline(ss, split, ':'))
      parts.push_back(split);

    // get chromosome and remaining string
    std::chrom = parts[0];
    std::start_end = parts[1];
    
    // reset and store region locations
    parts.clear();
    ss.str(start_end);
   
    // read in rest of string that shows region location 
    while(std::getline(ss, split, '-'))
      parts.push_back(split);

    // set first position to read from and set ending position
    start = stoi(parts[0]);
    reg_end = stoi(parts[1]);

    // end of first bin size
    end = stoi(start) + binsize; // - 1? TODO is it inclusive or exclusive?

    // check for size of region and ensure bin is inside
    if (end > reg_end)
      end = reg_end;
  }

  // first bin
  currentBin = new GenomeBin(chrom, start, end);
}


/*
  Inputs: none

  Outputs: 
  - bool: return true if there are still more bins left, else false

  Set currentBin to be the next bin
 */
bool BinGenerator::GotoNextBin() {
  // TODO Is the "end" parameter of GenomeBin inclusive?
  if (currentBin->end != reg_end)
  {
    // update the start and end for the next bin
    int32_t start = currentBin->start + binsize;
    int32_t end = currentBin->end + binsize;
    std::string chrom = currentBin->chrom;
    delete currentBin;

    // TODO Implement the default case. What is the true end?
    // How do we know if we've reached the end of the last chromosome?

    // Check if the new bin is outside the end region
    if (end <= reg_end)
      currentBin = new GenomeBin(chrom, start, end);
    else
      currentBin = new GenomeBin(chrom, start, reg_end);

    return true;
  }
  else
    return false;
}

/*
  Inputs: none

  Outputs:
  - GenomeBin: the current genome bin that was set in GotoNextBin
 */
const GenomeBin BinGenerator::GetCurrentBin() {
  return *currentBin;
}

BinGenerator::~BinGenerator() {
  delete currentBin;
}
