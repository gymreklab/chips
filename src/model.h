#ifndef SRC_MODEL_H__
#define SRC_MODEL_H__

#include <string>
#include "src/options.h"

class ChIPModel {
 public:
  // Constructors/destructors
  ChIPModel();
  virtual ~ChIPModel();

  // Write output
  void WriteModel(const std::string& output_json_file);
  void PrintModel();

  // Setters/getters
  void SetFrag(float k, float theta) {frag_param_k=k; frag_param_theta=theta;}
  void GetFrag(float* k, float*theta) {*k=frag_param_k; *theta=frag_param_theta;}
  void SetF(float f) {f_pulldown=f;}
  void GetF(float* f) {*f=f_pulldown;}
  void SetS(float s) {s_pulldown=s;}
  void GetS(float* s) {*s=s_pulldown;}
  void SetPCR(float pcr) {pcr_rate=pcr;}
  void GetPCR(float* pcr) {*pcr=pcr_rate;}
  void UpdateParams(const Options& options);
  void UpdateOptions(Options& options);
  void ReadFromJSON(const std::string& model_file);

 private:
  float frag_param_k;
  float frag_param_theta;
  float f_pulldown;
  float s_pulldown;
  float pcr_rate;
};

#endif  // SRC_MODEL_H__
