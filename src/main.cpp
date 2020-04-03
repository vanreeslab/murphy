#include <iostream>

#include <mpi.h>
#include "grid.hpp"
#include "field.hpp"

int main(int argc, char** argv){
    MPI_Init(&argc,&argv);

    bool periodic[3] = {false, true, false};
    int  l[3]        = {1, 2, 3};

    Grid*  grid     = new Grid(2, periodic, l, MPI_COMM_WORLD, NULL);
    Field* gaussian = new Field("gaussian", 3);

    grid->AddField(gaussian);

    delete(gaussian);
    delete(grid);


    MPI_Finalize();
}