#include "src/bingenerator.h"
#include "src/ref_genome.h"
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
 
  // default case
  if (options.region.empty())
  {
    // get the chroms and lengths from fasta file
    RefGenome ref (options.reffa);

    if (!ref.GetChroms(&chroms))
      PrintMessageDieOnError("Could not gather chromosomes from "
                                 + options.reffa, M_ERROR);
    if (!ref.GetLengths(&chromLengths))
      PrintMessageDieOnError("Could not gather chromosome lengths from "
                                 + options.reffa, M_ERROR);

    // index for next chromosome
    nextChrom = 1;
      
    // get last chromosome and the end of the first region
    endChrom = chroms[chroms.size() - 1];
    regEnd = chromLengths[chroms[0]];

    // set up parameters for GenomeBin
    chrom = chroms[0];
    start = 1;
    end = binsize;

    if (binsize > regEnd)
      end = regEnd;
  }

  // specified region
  else
  {
    vector<string> parts;            // vector to store split strings
    stringstream ss(options.region); // read in region
    string split;                    // used to store split strings

    // get chrom and locations from string
    while(getline(ss, split, ':'))
      parts.push_back(split);

    if (parts.size() != 2)
      PrintMessageDieOnError("Improper region input format should be chrom:start-end", M_ERROR);

    // get chromosome and remaining string
    chrom = parts[0];
    string start_end = parts[1];
    
    // set last chromosome to stop binning
    endChrom = chrom;
 
    // reset and store region locations
    parts.clear();
    ss.str(start_end);
    ss.clear();
   
    // read in rest of string that shows region location 
    while(getline(ss, split, '-'))
      parts.push_back(split);

    if (parts.size() != 2)
      PrintMessageDieOnError("Improper region input format should be chrom:start-end", M_ERROR);

    // set first position to read from and set ending position
    start = stoi(parts[0]);
    regEnd = stoi(parts[1]);

    // end of first bin size
    end = start + binsize - 1; // NOTE: It is inclusive meaning the end position is included.

    // check for size of region and ensure bin is inside
    if (end > regEnd)
      end = regEnd;
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

  // Go to the next bin
  if (currentBin->end != regEnd || currentBin->chrom != endChrom)
  {
    // update the start and end
    string chrom = currentBin->chrom;
    int start = currentBin->start + binsize;
    int end = currentBin->end + binsize;
    
    // update chromosome
    if (currentBin->end == regEnd && currentBin->chrom != endChrom)
    {
      chrom = chroms[nextChrom];
      start = 1;
      end = binsize;
      regEnd = chromLengths[chrom];
      nextChrom++;
    }

    delete currentBin;

    // Check if the new bin is outside the end region
    if (end <= regEnd)
      currentBin = new GenomeBin(chrom, start, end);
    else
      currentBin = new GenomeBin(chrom, start, regEnd);

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
