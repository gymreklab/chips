#include "src/model.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <json.hpp>

using namespace std;
using json = nlohmann::json;

ChIPModel::ChIPModel() {
  frag_param_k = -1;
  frag_param_theta = -1;
  f_pulldown = -1;
  s_pulldown = -1;
  pcr_rate = -1;
}

void ChIPModel::WriteModel(const std::string& output_json_file) {
  // Remove file if already there
  remove(output_json_file.c_str());

  // Set up writer
  json mwriter;

  // Fill in values
  mwriter["frag"]["k"] = frag_param_k;
  mwriter["frag"]["theta"] = frag_param_theta;
  mwriter["pulldown"]["f"] = f_pulldown;
  mwriter["pulldown"]["s"] = s_pulldown;
  mwriter["pcr_rate"] = pcr_rate;

  // Write to file
  ofstream outfile;
  outfile.open(output_json_file, ios_base::app);
  outfile << std::setw(4) << mwriter << std::endl;
  outfile.close();
}

void ChIPModel::PrintModel() {
  std::cerr << "K: " << frag_param_k << " Theta: " << frag_param_theta << endl;
  std::cerr << "f: " << f_pulldown << endl;
  std::cerr << "s: " << s_pulldown << endl;
  std::cerr << "pcr rate: " << pcr_rate << endl;
}

/*
  For any param that was set to -1, use the default
 */
void ChIPModel::UpdateParams(const Options& options) {
  if (frag_param_k == -1) {
    frag_param_k = options.gamma_k;
  }
  if (frag_param_theta == -1) {
    frag_param_theta = options.gamma_theta;
  }
  if (f_pulldown == -1) {
    f_pulldown = options.ratio_f;
  }
  if (s_pulldown == -1) {
    s_pulldown = options.ratio_s;
  }
  if (pcr_rate == -1) {
    pcr_rate = options.pcr_rate;
  }
}

/*
  Fill in options with model parameters
 */
void ChIPModel::UpdateOptions(Options& options) {
  options.pcr_rate = pcr_rate;
  options.gamma_k = frag_param_k;
  options.gamma_theta = frag_param_theta;
  options.ratio_f = f_pulldown;
  options.ratio_s = s_pulldown;
}

/*
  Update model params from JSON file.
  Don't overwrite anything already set (not -1)
 */
void ChIPModel::ReadFromJSON(const std::string& model_file) {
  json mjson;
  std::ifstream mreader(model_file);
  mreader >> mjson;
  if (mjson.find("pcr_rate") != mjson.end() &&
      pcr_rate == -1) {
    pcr_rate = mjson["pcr_rate"].get<double>();
  }
  if (mjson.find("pulldown") != mjson.end()) {
    if (mjson["pulldown"].find("f") != mjson["pulldown"].end() &&
	f_pulldown == -1) {
      f_pulldown = mjson["pulldown"]["f"].get<double>();
    }
    if (mjson["pulldown"].find("s") != mjson["pulldown"].end() &&
	s_pulldown == -1) {
      s_pulldown = mjson["pulldown"]["s"].get<double>();
    }
  }
  if (mjson.find("frag") != mjson.end()) {
    if (mjson["frag"].find("k") != mjson["frag"].end() &&
	frag_param_k == -1) {
      frag_param_k = mjson["frag"]["k"].get<double>();
    }
    if (mjson["frag"].find("s") != mjson["frag"].end() &&
	frag_param_theta == -1) {
      frag_param_theta = mjson["frag"]["s"].get<double>();
    }
  }
}

ChIPModel::~ChIPModel() {}
