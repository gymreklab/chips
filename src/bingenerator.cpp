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
  binsize = options.binsize;
  string chrom;
  int start;
  int end;
  
  if (options.region.empty() || options.region == "genome-wide")
  {
    // default is genome wide TODO figure out if its null or "genome-wide" as default
    // TODO How to handle traversing through whole genome?
    //      - How to check if we've reached the end of the first second third etc chr. 
  }
  else
  {
    vector<string> parts;        // vector to store split strings
    stringstream ss(options.region); // read in region
    string split;                     // used to store split strings

    // get chrom and locations from string
    while(getline(ss, split, ':'))
      parts.push_back(split);

    // get chromosome and remaining string
    chrom = parts[0];
    string start_end = parts[1];
    
    // reset and store region locations
    parts.clear();
    ss.str(start_end);
    ss.clear();
   
    // read in rest of string that shows region location 
    while(getline(ss, split, '-'))
      parts.push_back(split);

    // set first position to read from and set ending position
    start = stoi(parts[0]);
    reg_end = stoi(parts[1]);

    // end of first bin size
    end = start + binsize; // - 1? TODO is it inclusive or exclusive?

    // check for size of region and ensure bin is inside
    if (end > reg_end)
      end = reg_end;
  }
    
  // first bin
  firstBin = true;
  currentBin = new GenomeBin(chrom, start, end);
}


/*
  Inputs: none

  Outputs: 
  - bool: return true if there are still more bins left, else false

  Set currentBin to be the next bin
 */
bool BinGenerator::GotoNextBin() {
  if (firstBin)
  {
    firstBin = false;
    return true;
  }

  // TODO Is the "end" parameter of GenomeBin inclusive?
  if (currentBin->end != reg_end)
  {
    // update the start and end for the next bin
    int start = currentBin->start + binsize;
    int end = currentBin->end + binsize;
    string chrom = currentBin->chrom;
    delete currentBin;

    cout << "Updated Start: " << start << " Updated end: " << end << endl;
    cout << "Region end: " << reg_end << endl;

    // TODO Implement the default case. What is the true end?
    // Start at first chromosome and go where?
    // How do we know if we've reached the end of the any chromosome?

    // Check if the new bin is outside the end region
    if (end <= reg_end)
    {
      cout << "End is less or equal to reg_end" << endl;
      currentBin = new GenomeBin(chrom, start, end);
    }
    else
    {
      cout << "End was larger than End of region using end of region" << endl;
      currentBin = new GenomeBin(chrom, start, reg_end);
    }
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
