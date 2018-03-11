#ifndef SRC_SEQUENCER_H__
#define SRC_SEQUENCER_H__

#include "src/fragment.h"
#include "src/options.h"

#include <vector>

using namespace std;

class Sequencer {
 public:
  Sequencer(const Options& options);
  virtual ~Sequencer();

  void Sequence(const vector<Fragment>& input_fragments);
};

#endif  // SRC_SEQUENCER_H__
