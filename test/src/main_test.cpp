#include <mpi.h>

#include "gtest/gtest.h"
#include "murphy.hpp"

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;
    sc_init(comm, 1, 1, NULL, SC_LP_ESSENTIAL);
    p4est_init(NULL, SC_LP_PRODUCTION);
    //-------------------------------------------------------------------------
    {
        int comm_size;
        MPI_Comm_size(MPI_COMM_WORLD,&comm_size);
        if(comm_size != 1){
            m_log("unable to run tests with more than 1 CPU (for the moment)");
            MPI_Abort(MPI_COMM_WORLD,MPI_ERR_ASSERT);
        }
        ::testing::InitGoogleTest(&argc, argv);
        return RUN_ALL_TESTS();
    }
    //-------------------------------------------------------------------------
    sc_finalize();
    MPI_Finalize();
}