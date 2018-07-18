#ifndef SRC_OPTIONS_H__
#define SRC_OPTIONS_H__

#include <vector>
#include <string>

#include <stdint.h>

class Options {
 public:
  Options();

  // Required simulation inputs
  std::string reffa;
  std::string regionsbed;
  std::string outprefix;

  // Required learn inputs
  std::string chipbam;

  virtual ~Options();
};

#endif  // SRC_OPTIONS_H__
