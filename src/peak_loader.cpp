/* converter
 *
 * convert peakfile into a standard format:
 *      chromID  start  length  count
 * and save into a vector
 *
 * Last modified by An Zheng
 */
#include "src/peak_loader.h"

const std::map<std::string, int> PeakLoader::peakfileTypeList = {{"homer", 0}, {"test", 1},{"bed", 2}};

PeakLoader::PeakLoader(const std::string _peakfile, const std::string _peakfileType,
                        const std::string _bamfile, const std::int32_t _count_colidx){
  if(_peakfileType == ""){
    std::cerr << "****** ERROR: Need to specify the type of the peak file ******" << std::endl;
    std::exit(1);
  }else if ((_bamfile == "") && (_count_colidx == -1)){
    std::cerr << "****** ERROR: Need to specify either parameter -c or parameter -b ******" << std::endl;
    std::exit(1);
  }else{
    peakfile = _peakfile;
    peakfileType = _peakfileType;
    bamfile = _bamfile;
    count_colidx = _count_colidx;

    // unify names of peak file types into lower cases
    std::transform(_peakfileType.begin(), _peakfileType.end(), peakfileType.begin(), ::tolower);
  }
}

bool PeakLoader::Load(std::vector<Fragment>& peaks, const std::string region){
  PeakReader peakreader(peakfile);
  switch (peakfileTypeList.at(peakfileType)){
    case 0:
      peakreader.HomerPeakReader(peaks, count_colidx, region);
      break;
    case 1:
      peakreader.TestPeakReader(peaks, count_colidx, region);
      break;
    case 2:
      peakreader.BedPeakReader(peaks, count_colidx, region);
      break;
    default:
      std::cerr << "An unexpected error happened in PeakLoader->Load()" << std::endl;
      return false;
      break;
  }
  if (bamfile != ""){
    peakreader.UpdateTagCount(peaks, bamfile, &total_genome_length, &total_tagcount, &tagcount_in_peaks, region);
  }
  return true;
}
