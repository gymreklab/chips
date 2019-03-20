#include "src/library_constructor.h"

LibraryConstructor::LibraryConstructor(const Options& options) {
  pcr_rate = options.pcr_rate;
}

// Apply PCR
void LibraryConstructor::Perform(const vector<Fragment>& input_fragments,
				 vector<Fragment>* output_fragments) {
  for (int frag_index=0; frag_index<input_fragments.size(); frag_index++){
    while (true) {
      output_fragments->push_back(input_fragments[frag_index]);
      if (rand()/double(RAND_MAX) < pcr_rate) break;
    }
  }
}


LibraryConstructor::~LibraryConstructor() {}
