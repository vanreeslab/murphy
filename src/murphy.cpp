#include <argp.h>
#include <mpi.h>

#include <string>

#include "clients/flow_abc.hpp"
#include "clients/navier_stokes.hpp"
#include "clients/simple_advection.hpp"
#include "core/macros.hpp"
#include "core/types.hpp"
#include "p8est.h"
#include "parser.hpp"

using std::string;

TestCase* MurphyInit(int argc, char** argv) {
    //-------------------------------------------------------------------------
    int provided;
    // set MPI_THREAD_SERIALIZED
    int requested = MPI_THREAD_SERIALIZED;
    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;

    // sc_init(comm, 1, 1, NULL, SC_LP_SILENT);
    sc_init(comm, 1, 1, NULL, SC_LP_INFO);
    p4est_init(NULL, SC_LP_SILENT);
    // p4est_init(NULL, SC_LP_INFO);

    // so dome checks for the aligment, the constants etc
    m_assert(M_GS >= 1, "1 is the min ghost point needed, because of the IO");
    m_assert(M_N >= M_GS, "we cannot have ghost points that span more than 1 block");
    m_assert((M_STRIDE * M_GS + M_GS) % (M_ALIGNMENT / sizeof(real_t)) == 0, "the first point has to be aligned");
    m_assert((M_STRIDE*M_STRIDE*M_STRIDE) < std::numeric_limits<bidx_t>::max(),"The whole block indexing must be contained in bidx_t");

    // parse arguments, just simply display the help, nothing more at that stage
    ParserArguments argument;
    ParseArgument(argc, argv, &argument);

    // init the testcase
    TestCase* testcase;
    if (argument.do_navier_stokes) {
        testcase = new NavierStokes();
        testcase->InitParam(&argument);
        return testcase;
    } 
    else if (argument.do_abc_flow){
        testcase = new FlowABC();
        testcase->InitParam(&argument);
        return testcase;
    }
    else if (argument.do_simple_adv){
        testcase = new SimpleAdvection();
        testcase->InitParam(&argument);
        return testcase;
    }
    else {
        // m_assert(false, "no testcase has been choosen");
        m_log("no test case chosen");
    }
    return nullptr;
    //-------------------------------------------------------------------------
}

void MurphyFinalize(TestCase* testcase) {
    //-------------------------------------------------------------------------
    if (testcase != nullptr) {
        delete (testcase);
    }
    m_log("bye bye MURPHY");

    sc_finalize();
    MPI_Finalize();
    //-------------------------------------------------------------------------
}