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
        Grid* grid = new Grid(argument.init_lvl_, argument.period_, argument.length_, MPI_COMM_WORLD, prof);

        // get an refined and adapted grid given the patch
        list<Patch> patch;
        real_t      origin[3] = {0.5, 0.0, 0.0};
        real_t      length[3] = {0.5, 1.0, 1.0};
        patch.push_back(Patch(origin, length, 2));
        grid->Adapt(&patch);
        grid->Adapt(&argument.patch_);
        // create a field
        Field* vort = new Field("vorticity", 3);
        // Field* psi  = new Field("psi", 3);
        // Field* res  = new Field("residual", 3);
        // // register the field to the grid
        grid->AddField(vort);
        // grid->AddField(psi);
        // grid->AddField(res);
        // // get some values for the polynomial
        // real_t     dir[3]  = {0.0, 1.0, 0.0};
        // lid_t      deg[3]  = {0, 1, 0};
        // SetPolynom polynom = SetPolynom(deg, dir);
        // polynom(grid, vort);
        // real_t length[3] = {1.0* argument.length_[0],1.0* argument.length_[1],1.0* argument.length_[2]};
        // real_t freq[3] = {2.0,4.0,1.0};
        // SetSinus sinus = SetSinus(length,freq);
        // sinus(grid,vort);

        const real_t  center[3] = {0.5, 0.5, 0.5};
        const lda_t   normal    = 2;
        const real_t  sigma     = 0.2;
        const real_t  radius    = 0.25;
        SetVortexRing vr_init(normal, center, sigma, radius);
        vr_init(grid, vort);

        // // set an EVEN bc for everybody (everywhere and in X direction for each dimension)
        // vort->bctype(M_BC_ODD);
        vort->bctype(M_BC_EXTRAP_5);
        // psi->bctype(M_BC_ODD);
        // res->bctype(M_BC_ODD);

        IOH5 dump = IOH5("data");
        dump(grid, vort);

        // grid->Adapt(vort);

        // grid->GhostPull(vort);
        // dump.dump_ghost(true);
        // dump(grid, vort);

        // init the MG solver
        // Multigrid* poisson = new Multigrid(grid, 0, vort, psi, res);

        // poisson->Solve();

        // and destroy the grid and the field
        delete (vort);
        // delete (psi);
        // delete (res);
        // delete (poisson);
        delete (grid);
    }
    // display the profiler
    prof->Disp();
    delete (prof);
    m_log("leaving, bye bye murphy\n");
    //-------------------------------------------------------------------------
    murphy_finalize();
}