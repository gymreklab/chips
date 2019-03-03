#include <iostream>

#include "src/common.h"

using namespace std;

// Function declarations
int chipmunk_help(void);
int simulate_reads_main(int argc, char* argv[1]);
int learn_main(int argc, char* argv[1]);

int main(int argc, char* argv[]) {
  // make sure the user at least entered a sub_command
  if (argc < 2) return chipmunk_help();

  std::string sub_cmd = argv[1];

  if (sub_cmd == "simreads") {
    return simulate_reads_main(argc-1, argv+1);
  } else if (sub_cmd == "learn") {
    return learn_main(argc-1, argv+1);
  } else if (sub_cmd == "-h" || sub_cmd == "--help" ||
	   sub_cmd == "-help") {
    return chipmunk_help();
  } else if (sub_cmd == "-version" || sub_cmd == "--version") {
    cout << "chipmunk " << _GIT_VERSION << endl;
  } else {
    cerr << "error: unrecognized command: " << argv[1] << endl << endl;
    return 1;
  }
  return 0;
}

int chipmunk_help(void) {
  cout << PROGRAM_NAME << ": Simulator for ChIP-seq and other -seq experiments.\n";
  cout << "usage:   chipmunk <subcommand> [options]" << endl << endl;
  cout << "The chipmunk sub-commands include:"<<endl;
  cout << "     simreads      " << "Simulate ChIP-seq reads given a set of intervals.\n";
  cout << "     learn         " << "Learn model from real ChIP data.\n";
  cout  << endl;
  cout  << "[ General help ]" << endl;
  cout  << "    --help        "  << "Print this help menu.\n";
  cout  << "    --version     "  << "What version are you using?.\n";
  return 0;
}
