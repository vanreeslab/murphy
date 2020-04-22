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
        omp_set_num_threads(1);

        bool periodic[3] = {false, false, false};
        int  l[3]        = {2, 2, 2};
        // create a grid
        Grid* grid = new Grid(0, periodic, l, MPI_COMM_WORLD, NULL);
        // create a field
        Field* vort = new Field("vorticity", 3);
        // register the field to the grid
        grid->AddField(vort);

        // set a Gaussian
        real_t center[3] = {l[0] * 0.25, l[1] * 0.5, l[2] * 0.5};
        real_t sigma[3]  = {l[0] * 0.1, l[1] * 0.0, l[2] * 0.0};
        real_t freq[3]   = {0.0, 0.5, 0.5};
        real_t length[3] = {(real_t)l[0], (real_t)l[1], (real_t)l[2]};

        SetExpoCosinus expocos = SetExpoCosinus(center, sigma, length, freq);
        expocos(grid, vort);
        // set an EVEN bc for everybody (everywhere and in X direction for each dimension)
        vort->bctype(M_BC_EVEN);
        // for(int id=0; id<3; id++){
        //     vort->bctype(M_BC_EXTRAP,id,0);
        //     vort->bctype(M_BC_EXTRAP,id,1);
        // }

        // get an refined and adapted grid
        grid->Adapt(vort);

        // create a dumper and dump
        IOH5 mydump = IOH5("data");
        mydump(grid, vort,"vort_fine");

        grid->GhostPull(vort);

        mydump.dump_ghost(true);
        mydump(grid, vort,"vort_fine");

        // destroy the grid and the field
        delete (vort);
        delete (grid);
    }

    m_log("leaving, bye bye murphy\n");
    //-------------------------------------------------------------------------
    sc_finalize();
    MPI_Finalize();
}