#include <mpi.h>

#include <iostream>
#include <string>

#include "murphy.hpp"
#include "testcase.hpp"

using std::string;
using std::to_string;

int main(int argc, char** argv) {
    //-------------------------------------------------------------------------
    TestCase* testcase = MurphyInit(argc, argv);

    if (testcase != nullptr) {
        testcase->Run();
    }

    MurphyFinalize(testcase);
    //-------------------------------------------------------------------------
}
