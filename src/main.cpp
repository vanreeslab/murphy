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
#include "parser.hpp"

int main(int argc, char** argv) {
    murphy_init(argc,argv);
    //-------------------------------------------------------------------------
    {
        parse_arg_t argument;
        ParseArgument(argc, argv, &argument);
        // create a grid
        Prof* prof = new Prof("MURPHY");
        Grid* grid = new Grid(1, argument.period_,argument.length_, MPI_COMM_WORLD, prof);
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
        // IOH5 dump = IOH5("data");

        // get an refined and adapted grid
        // grid->Adapt(vort);
        grid->Adapt(&argument.patch_);

        // dump(grid, vort);

        LaplacianCross<5> lapla = LaplacianCross<5>(grid);
        lapla(vort,diff);

        // dump(grid,diff);

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
    murphy_finalize();
}