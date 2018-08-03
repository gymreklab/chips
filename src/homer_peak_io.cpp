// Taken from https://github.com/tfwillems/HipSTR/blob/master/src/bam_io.cpp


#include "homer_peak_io.h"

HomerPeakReader::HomerPeakReader(const std::string& peakfile, std::vector<seq_loc>& peaks){
  std::string line;
  std::ifstream infile(peakfile.c_str());
  while (std::getline(infile, line)){
    if (line.find('#') != 0){
      std::string peakID;
      std::string chr;
      uint32_t start;
      uint32_t end;

      std::stringstream linestream(line);
      linestream >> peakID >> chr >> start >> end;
      
      seq_loc location = {chr, start, end};
      peaks.push_back(location);
    }
  }
}
/* END OF FILE */
