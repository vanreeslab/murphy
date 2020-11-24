#ifndef SRC_MURPHY_HPP_
#define SRC_MURPHY_HPP_

#include <p8est.h>

#include <limits>
#include <map>
#include <string>

#include "defs.hpp"
#include "testcase.hpp"

TestCase* MurphyInit(int argc, char* argv[]);
void      MurphyFinalize(TestCase* testcase);

#endif  // SRC_MURPHY_HPP_
