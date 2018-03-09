#include <cstring>
#include <iostream>
#include <stdlib.h>

using namespace std;

// define our program name
#define PROGRAM_NAME "asimon simreads"

// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

// Function declarations
void simulate_reads_help(void);

int simulate_reads_main(int argc, char* argv[]) {
  bool showHelp = false;

  // check to see if we should print out some help
  if(argc <= 1) showHelp = true;
  for(int i = 1; i < argc; i++) {
    int parameterLength = (int)strlen(argv[i]);

    if((PARAMETER_CHECK("-h", 2, parameterLength)) ||
       (PARAMETER_CHECK("--help", 5, parameterLength))) {
      showHelp = true;
    }
  }
  if (showHelp) {simulate_reads_help();}

  // TODO parse other params

  if (!showHelp) {
    // TODO do something!
  } else {
    simulate_reads_help();
    return 0;
  }
}

void simulate_reads_help(void) {
  cerr << "\nTool:    asimon simreads" << endl;
  cerr << "Version: " << _GIT_VERSION << "\n";    
  cerr << "Summary: Simulate ChIP-seq reads for a set of peaks." << endl << endl;
  cerr << "Usage:   " << PROGRAM_NAME << " TODO " << endl << endl;
  // end the program here
  exit(1);
}
