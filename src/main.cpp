#include <mpi.h>

#include <iostream>

#include "field.hpp"
#include "grid.hpp"
#include "ioh5.hpp"
#include "setvalues.hpp"
#include "ghost.hpp"

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;
    sc_init(comm, 1, 1, NULL, SC_LP_ESSENTIAL);
    p4est_init(NULL, SC_LP_PRODUCTION);
    //-------------------------------------------------------------------------
    {
        bool periodic[3] = {false, true, false};
        int  l[3]        = {1, 2, 3};

        // create a grid
        Grid*  grid = new Grid(1, periodic, l, MPI_COMM_WORLD, NULL);
        // create a field
        Field* vort = new Field("vorticity", 3);
        grid->AddField(vort);
        // set a Gaussian
        real_t      center[3] = {l[0] * 0.5, l[1] * 0.5, l[2] * 0.5};
        SetGaussian gaussian  = SetGaussian(0.1, center);
        gaussian(grid, vort);
        // create a dumper and dump
        // IOH5 mydump = IOH5("data");
        // mydump(grid, vort);

        // create a ghost 
        Ghost* ghost = new Ghost(grid);

    
        // destroy the grid and the field
        delete(ghost);
        delete (vort);
        delete (grid);
    }
    //-------------------------------------------------------------------------
    sc_finalize();
    MPI_Finalize();
}