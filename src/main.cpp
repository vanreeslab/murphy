#include <mpi.h>

#include <iostream>

#include "field.hpp"
#include "grid.hpp"
#include "ioh5.hpp"
#include "setvalues.hpp"
#include "laplacian.hpp"
#include "prof.hpp"
#include "parser.hpp"

int main(int argc, char** argv) {
    murphy_init(argc, argv);
    //-------------------------------------------------------------------------
    {
        parse_arg_t argument;
        ParseArgument(argc, argv, &argument);
        // create a grid
        Prof* prof = new Prof("MURPHY");
        Grid* grid = new Grid(0, argument.period_, argument.length_, MPI_COMM_WORLD, prof);
        // get an refined and adapted grid given the patch
        grid->Adapt(&argument.patch_);
        // create a field
        Field* vort = new Field("vorticity", 3);
        Field* diff = new Field("diffusion", 3);
        // register the field to the grid
        grid->AddField(vort);
        grid->AddField(diff);
        // get some values for the polynomial
        real_t     dir[3]  = {1.0, 1.0, 1.0};
        lid_t      deg[3]  = {2, 2, 2};
        SetPolynom polynom = SetPolynom(deg, dir);
        polynom(grid, vort);
        // set an EVEN bc for everybody (everywhere and in X direction for each dimension)
        vort->bctype(M_BC_EXTRAP_5);
        // take the derivtive and store it into the diff field
        LaplacianCross<5> lapla = LaplacianCross<5>(grid);
        lapla(vort, diff);

        // display the profiler
        prof->Disp();
        // and destroy the grid and the field
        delete (vort);
        delete (diff);
        delete (grid);
        delete (prof);
    }

    m_log("leaving, bye bye murphy\n");
    //-------------------------------------------------------------------------
    murphy_finalize();
}