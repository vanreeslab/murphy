#include <argp.h>
#include <mpi.h>

#include <string>

#include "clients/convergence_weno.hpp"
#include "clients/debug_lifting.hpp"
#include "clients/epsilon_test.hpp"
#include "clients/flow_abc.hpp"
#include "clients/navier_stokes.hpp"
#include "clients/simple_advection.hpp"
#include "clients/twolevel_convweno.hpp"
#include "clients/weak_scalability.hpp"
#include "core/macros.hpp"
#include "core/types.hpp"
#include "p8est.h"
#include "parser.hpp"

using std::string;

/**
 * @brief writes the file murphy.info used for tracking of the results, bookkeeping etc
 */
void WriteInfo(int argc, char** argv) {
    rank_t rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        string filename(argv[0]);
        // remove the "./" that might be in it
        filename.erase(std::remove(filename.begin(), filename.end(), '.'), filename.end());
        filename.erase(std::remove(filename.begin(), filename.end(), '/'), filename.end());
        filename += ".info";
        FILE* file = fopen(filename.c_str(), "w+");
        fprintf(file, "MURPHY - (c) MIT\n");
        fprintf(file, "- commit: %s\n", M_GIT_COMMIT);
        fprintf(file, "- defines:\n");
        fprintf(file, "\tM_N = %d\n", M_N);
        fprintf(file, "\tM_GS = %d\n", M_GS);
        fprintf(file, "\tM_WAVELET_N = %d\n", M_WAVELET_N);
        fprintf(file, "\tM_WAVELET_NT = %d\n", M_WAVELET_NT);
        fprintf(file, "- argument list:\n");
        for (int i = 0; i < argc; ++i) {
            fprintf(file, "\t%s\n", argv[i]);
        }
        fclose(file);
    }
}

TestCase* MurphyInit(int argc, char** argv) {
    //-------------------------------------------------------------------------
    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;
    sc_init(comm, 0, 0, NULL, SC_LP_SILENT);
    p4est_init(NULL, SC_LP_SILENT);

#ifdef M_MPI_AGGRESSIVE
    // set some properties on the COMM_WORLD to be more aggressive
    MPI_Info info;
    MPI_Info_create(&info);
    MPI_Info_set(info, "mpi_assert_exact_length", "true");
    MPI_Info_set(info, "mpi_assert_allow_overtaking", "true");
    MPI_Comm_set_info(MPI_COMM_WORLD, info);
    MPI_Info_free(&info);
#endif

    // so dome checks for the aligment, the constants etc
    m_assert((M_N % 2) == 0, "the number of points must be odd");
    m_assert(M_GS >= 1, "1 is the min ghost point needed, because of the IO");
    m_assert(M_N >= M_GS, "we cannot have ghost points that span more than 1 block");
    // m_assert((M_STRIDE * M_GS + M_GS) % (M_ALIGNMENT / sizeof(real_t)) == 0, "the first point has to be aligned");

    // write the info file
    WriteInfo(argc, argv);

    // parse arguments, just simply display the help, nothing more at that stage
    ParserArguments argument;
    ParseArgument(argc, argv, &argument);

    // init the testcase
    TestCase* testcase;
    if (argument.do_navier_stokes) {
        testcase = new NavierStokes();
        testcase->InitParam(&argument);
        return testcase;
    } else if (argument.do_abc_flow) {
        testcase = new FlowABC();
        testcase->InitParam(&argument);
        return testcase;
    } else if (argument.do_simple_adv) {
        testcase = new SimpleAdvection();
        testcase->InitParam(&argument);
        return testcase;
    } else if (argument.do_epsilon_test) {
        testcase = new EpsilonTest();
        testcase->InitParam(&argument);
        return testcase;
    } else if (argument.do_debug_lifting) {
        testcase = new DebugLifting();
        testcase->InitParam(&argument);
        return testcase;
    } else if (argument.do_conv_weno) {
        testcase = new ConvergenceWeno();
        testcase->InitParam(&argument);
        return testcase;
    } else if (argument.do_2lvl_weno) {
        testcase = new TwoLevelConvWeno();
        testcase->InitParam(&argument);
        return testcase;
    } else if (argument.do_weak_scal) {
        testcase = new WeakScalability();
        testcase->InitParam(&argument);
        return testcase;
    } else {
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
