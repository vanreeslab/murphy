#include <mpi.h>

#include <iostream>
#include <string>

#include "advection_diffusion.hpp"
#include "field.hpp"
#include "grid.hpp"
#include "ioh5.hpp"
#include "laplacian.hpp"
#include "multigrid.hpp"
#include "murphy.hpp"
#include "parser.hpp"
#include "prof.hpp"
#include "setvalues.hpp"

using std::string;
using std::to_string;

int main(int argc, char** argv) {
    //-------------------------------------------------------------------------
    TestCase* testcase = MurphyInit(argc, argv);

    testcase->Run();

    MurphyFinalize(testcase);
    //-------------------------------------------------------------------------
}
