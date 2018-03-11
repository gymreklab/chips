#ifndef SRC_SHEARER_H__
#define SRC_SHEARER_H__

#include "src/fragment.h"
#include "src/options.h"

#include <vector>

using namespace std;

class Shearer {
 public:
  Shearer(const Options& _options);
  void Shear(vector<Fragment>* sheared_fragments);
  virtual ~Shearer();
};

#endif  // SRC_SHEARER_H__
