#ifndef SRC_PULLDOWN_H__
#define SRC_PULLDOWN_H__

#include "src/fragment.h"
#include "src/options.h"

#include <vector>

using namespace std;

class Pulldown {
 public:
  Pulldown(const Options& options);
  virtual ~Pulldown();

  void Perform(const vector<Fragment>& input_fragments,
	       vector<Fragment>* output_fragments);
};

#endif  // SRC_PULLDOWN_H__
