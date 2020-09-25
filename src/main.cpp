#include <mpi.h>

#include <iostream>
#include <string>

#include "field.hpp"
#include "grid.hpp"
#include "ioh5.hpp"
#include "setvalues.hpp"
#include "laplacian.hpp"
#include "prof.hpp"
#include "parser.hpp"
#include "multigrid.hpp"

using std::string;
using std::to_string;

int main(int argc, char** argv) {
    murphy_init(argc, argv);
    // get the argument lists
    parse_arg_t argument;
    ParseArgument(argc, argv, &argument);

    int comm_size;
    MPI_Comm_size(MPI_COMM_WORLD,&comm_size);
    string prof_name = string("MURPHY_") + to_string(comm_size) + string("ranks_") + to_string(omp_get_max_threads()) + string("threads");
    Prof*  prof      = new Prof(prof_name);

    //-------------------------------------------------------------------------
    for (int i = 0; i < argument.n_repeat_; i++) {
        // create a grid
        m_log("init level? %d", argument.init_lvl_);
        m_log("number of roots? %d %d %d", argument.length_[0], argument.length_[1], argument.length_[2]);
        m_log("periodic? %d %d %d", argument.period_[0], argument.period_[1], argument.period_[2]);
        
        // create the grid, at this point, no memory is created
        Grid grid(argument.init_lvl_, argument.period_, argument.length_, MPI_COMM_WORLD, prof);

        // // get an refined and adapted grid given the patch
        // list<Patch> patch;
        // real_t      origin[3] = {0.5, 0.0, 0.0};
        // real_t      length[3] = {0.5, 1.0, 1.0};
        // // patch.push_back(Patch(origin, length, 2));
        // grid->Adapt(&patch);
        // grid->Adapt(&argument.patch_);

        // create a field
        Field vort("vorticity", 3);
        grid.AddField(&vort);

        const real_t  center[3] = {argument.length_[0] / 2.0, argument.length_[1] / 2.0, argument.length_[2] / 2.0};
        const lda_t   normal    = 2;
        const real_t  sigma     = 0.05;
        const real_t  radius    = 0.25;
        SetVortexRing vr_init(normal, center, sigma, radius, grid.NGhostFront(), grid.NGhostBack());

        // set the BC for kiding
        vort.bctype(M_BC_EXTRAP_3);

        // adapt the mesh
        grid.SetTol(1e-1, 1e-2);
        grid.AdaptInitialCondition(&vort,&vr_init);

        // create the IO
        IOH5 dump("data");
        dump(&grid, &vort);
        grid.GhostPull(&vort);
        // dump.dump_ghost(true);
        // dump(&grid, &vort);

    }
    // display the profiler
    prof->Disp();
    delete (prof);
    m_log("leaving, bye bye murphy\n");
    //-------------------------------------------------------------------------
    murphy_finalize();
}
