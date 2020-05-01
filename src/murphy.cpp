#include "murphy.hpp"

#include <mpi.h>

#include "p8est.h"

void murphy_init(int argc, char** argv) {
    m_begin;
    //-------------------------------------------------------------------------
    int provided;
    // set MPI_THREAD_FUNNELED or MPI_THREAD_SERIALIZED
    int requested = MPI_THREAD_FUNNELED;
    MPI_Init_thread(&argc, &argv, requested, &provided);
    if (provided != requested) {
        printf("The MPI-provided thread behavior does not match\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    MPI_Comm comm = MPI_COMM_WORLD;
    sc_init(comm, 1, 1, NULL, SC_LP_SILENT);
    p4est_init(NULL, SC_LP_SILENT);
    // so dome checks for the aligment, the constants etc
    m_assert((M_STRIDE * M_GS + M_GS) % (M_ALIGNMENT / sizeof(real_t)) == 0, "the first point has to be aligned");
    //-------------------------------------------------------------------------
    m_end;
}

void murphy_finalize() {
    sc_finalize();
    MPI_Finalize();
}