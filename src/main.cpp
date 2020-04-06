#include <iostream>

#include <mpi.h>
#include "grid.hpp"
#include "field.hpp"
#include "setvalues.hpp"

int main(int argc, char** argv){
    MPI_Init(&argc,&argv);

    bool periodic[3] = {false, true, false};
    int  l[3]        = {1, 2, 3};

    Grid*  grid     = new Grid(2, periodic, l, MPI_COMM_WORLD, NULL);
    Field* vort = new Field("vorticity", 3);
    grid->AddField(vort);

    // set a Gaussian
    real_t center[3] = {l[0]*0.5,l[1]*0.5,l[2]*0.5};
    SetGaussian gaussian = SetGaussian(0.1,center);

    // DoOp<OperatorF*>(CallOp,grid,vort,&gaussian);
    gaussian.DoOp(grid,vort);


    delete(vort);
    delete(grid);


    MPI_Finalize();
}