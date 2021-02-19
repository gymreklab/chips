#include <err.h>
#include <stdlib.h>

#include <iostream>
#include <sstream>

#include "common.h"
#include "chipsConfig.h"

using namespace std;

void PrintMessageDieOnError(const string& msg, MSGTYPE msgtype) {
  string typestring = "";
  switch (msgtype) {
  case M_ERROR:
    typestring = "ERROR: ";
    break;
  case M_WARNING:
    typestring = "WARNING: ";
    break;
  case M_PROGRESS:
    typestring = "ProgressMeter: ";
    break;
  case M_DEBUG:
    typestring = "DEBUG: ";
    break;
  default:
    errx(1, "Invalid message type. This should never happen");
  }
  stringstream ss;
  ss << "[chips-" << chips_VERSION_MAJOR << "." << chips_VERSION_MINOR << "]"
     << typestring << msg << endl;
  cerr << ss.str();

  if (msgtype == M_ERROR) {
    exit(1);
  }
}

void RegionParser(const std::string region,
		  std::string& chromID, std::int32_t& start, std::int32_t& end){
  std::stringstream ss(region);

  // parse chromID and the rest
  std::string split;
  std::vector<std::string> parts;
  while(std::getline(ss, split, ':')) parts.push_back(split);
  if (parts.size() != 2) PrintMessageDieOnError("Improper region input format should be chrom:start-end", M_ERROR);
  chromID = parts[0];
  std::string start_end = parts[1];

  // reset
  parts.clear();
  ss.str(start_end);
  ss.clear();

  // parse start and end
  while(getline(ss, split, '-')) parts.push_back(split);
  if (parts.size() != 2) PrintMessageDieOnError("Improper region input format should be chrom:start-end", M_ERROR);
  start = stoi(parts[0]);
  end = stoi(parts[1]);
}
