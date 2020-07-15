#include <err.h>
#include <stdlib.h>

#include <iostream>
#include <sstream>

#include "common.h"

using namespace std;

void PrintMessageDieOnError(const string& msg, MSGTYPE msgtype) {
  string typestring = "";
  switch (msgtype) {
  case M_ERROR:
    typestring = "ERROR: ";
    break;
  case M_WARNING:
    typestring = "WARNING: ";
    break;
  case M_PROGRESS:
    typestring = "ProgressMeter: ";
    break;
  case M_DEBUG:
    typestring = "DEBUG: ";
    break;
  default:
    errx(1, "Invalid message type. This should never happen");
  }
  stringstream ss;
  ss  << "[" << PROGRAM_NAME
      << "-" << _GIT_VERSION << "] " << typestring << msg << endl;
  cerr << ss.str();

  if (msgtype == M_ERROR) {
    exit(1);
  }
}
