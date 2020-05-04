#include <mpi.h>

#include "gtest/gtest.h"
#include "murphy.hpp"

int main(int argc, char** argv) {
    // we should first init Google Test and then MPI.
    // for (int ia = 0; ia < argc; ia++) {
    //     printf("command = %s\n", argv[ia]);
    // }
    // but I need the rank to init google test, so I do the opposite
    // init MPI and pe4st
    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;
    sc_init(comm, 1, 1, NULL, SC_LP_ESSENTIAL);
    p4est_init(NULL, SC_LP_PRODUCTION);
    //-------------------------------------------------------------------------
    int err;
    {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        for (int ia = 0; ia < argc; ia++) {
            int comp;
            comp = strcmp(argv[ia], "--gtest_output=xml:report_test.xml");
            if (comp == 0) {
                sprintf(argv[ia], "--gtest_output=xml:report_test_rank%d.xml", rank);
            }

            comp = strcmp(argv[ia], "--gtest_output=xml:report_valid.xml");
            if (comp == 0) {
                sprintf(argv[ia], "--gtest_output=xml:report_valid_rank%d.xml", rank);
            }
        }
        // init google and remove Google Test arguments
        ::testing::InitGoogleTest(&argc, argv);
        err = RUN_ALL_TESTS();
    }
    //-------------------------------------------------------------------------
    sc_finalize();
    MPI_Finalize();

    return err;
}