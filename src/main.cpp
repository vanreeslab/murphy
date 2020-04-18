#include <mpi.h>

#include <iostream>

#include "field.hpp"
#include "ghost.hpp"
#include "grid.hpp"
#include "ioh5.hpp"
#include "setvalues.hpp"
#include "wavelet.hpp"

int main(int argc, char** argv) {
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
    //-------------------------------------------------------------------------
    {
        bool periodic[3] = {true, true, true};
        int  l[3]        = {3, 3, 3};

        // create a grid
        Grid* grid = new Grid(0, periodic, l, MPI_COMM_WORLD, NULL);
        // create a field
        Field* vort = new Field("vorticity", 3);
        grid->AddField(vort);
        // set a Gaussian
        real_t      center[3] = {l[0] * 0.5, l[1] * 0.5, l[2] * 0.5};
        SetGaussian gaussian  = SetGaussian(0.2, center);
        gaussian(grid, vort);

        // real_t alpha[3] = {0.0, 0.0, 1.0};
        // // SetAbs setabs   = SetAbs(alpha, center);
        // // setabs(grid, vort);
        // SetJump setjump   = SetJump(alpha, center);
        // setjump(grid, vort);

        // grid->GhostPull(vort);

        // create a dumper and dump
        // IOH5 mydump2 = IOH5("data");
        // mydump2.dump_ghost(true);
        // mydump2(grid, vort);
        // grid->Refine(1);

        // create a dumper and dump
        IOH5 mydump = IOH5("data");
        mydump(grid, vort,"vort_fine");

        grid->Adapt(vort);

        grid->GhostPull(vort);

        mydump(grid, vort,"vort_adapt");
        mydump.dump_ghost(true);
        mydump(grid, vort,"vort_adapt");


        // destroy the grid and the field
        delete (vort);
        delete (grid);
    }

    m_log("leaving, bye bye murphy\n");
    //-------------------------------------------------------------------------
    sc_finalize();
    MPI_Finalize();
}