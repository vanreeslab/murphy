#include <mpi.h>

#include <iostream>

#include "field.hpp"
#include "ghost.hpp"
#include "grid.hpp"
#include "ioh5.hpp"
#include "setvalues.hpp"
#include "wavelet.hpp"
#include "laplacian.hpp"
#include "prof.hpp"

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
        bool periodic[3] = {false, false, false};
        int  l[3]        = {3, 2, 1};
        // create a grid

        Prof* prof = new Prof("MURPHY");
        Grid* grid = new Grid(0, periodic, l, MPI_COMM_WORLD, prof);
        // create a field
        Field* vort = new Field("vorticity", 3);
        Field* diff = new Field("diffusion", 3);
        // register the field to the grid
        grid->AddField(vort);
        grid->AddField(diff);

        real_t     dir[3]  = {1.0, 1.0, 1.0};
        lid_t      deg[3]  = {2, 2, 2};
        SetPolynom polynom = SetPolynom(deg, dir);
        polynom(grid, vort);
        // set an EVEN bc for everybody (everywhere and in X direction for each dimension)
        vort->bctype(M_BC_EVEN);
        // set the extrapolation in the X direction [0,1]
        for(int id=0; id<3; id++){
            vort->bctype(M_BC_EXTRAP_3,id,0);
            vort->bctype(M_BC_EXTRAP_3,id,1);
        }
        //  // create a dumper and dump
        IOH5 dump = IOH5("data");
        // mydump(grid, vort,"vort");
        

        // get an refined and adapted grid
        grid->Adapt(vort);
        grid->Adapt(vort);

        dump(grid, vort);

        // grid->GhostPull(vort);

        // create a dumper and dump
        // dump.dump_ghost(true);
        // dump(grid, vort,);

        LaplacianCross<5> lapla = LaplacianCross<5>(grid);
        lapla(vort,diff);

        dump(grid,diff);

        // grid->GhostPull(vort);

        // mydump.dump_ghost(true);
        // mydump(grid, vort,"vort_fine");

        prof->Disp();

        // destroy the grid and the field
        delete (vort);
        delete (grid);
        delete (prof);
    }

    m_log("leaving, bye bye murphy\n");
    //-------------------------------------------------------------------------
    sc_finalize();
    MPI_Finalize();
}