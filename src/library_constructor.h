#ifndef SRC_LIBRARY_CONSTRUCTOR_H__
#define SRC_LIBRARY_CONSTRUCTOR_H__

#include "src/fragment.h"
#include "src/options.h"

#include <vector>

using namespace std;

class LibraryConstructor {
 public:
  LibraryConstructor(const Options& options);
  virtual ~LibraryConstructor();

  void Perform(const vector<Fragment>& input_fragments,
	       vector<Fragment>* output_fragments);
};

#endif  // SRC_LIBRARY_CONSTRUCTOR_H__
