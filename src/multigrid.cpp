#include "multigrid.hpp"

#include <mpi.h>
#include <limits>
#include <iostream>

#include "defs.hpp"
#include "gridcallback.hpp"
#include "laplacian.hpp"
#include "daxpy.hpp"
#include "gaussseidel.hpp"
#include "error.hpp"

#include "ioh5.hpp"

using std::numeric_limits;

Multigrid::Multigrid(Grid* grid, sid_t fft_level,Field* rhs, Field* sol,Field* res) {
    m_begin;
    //-------------------------------------------------------------------------
    // store the desired fft level and which field will be used as what
    fft_level_ = fft_level;
    m_assert(grid->IsAField(rhs), "the source field MUST exist in the grid");
    fields_nickname_["rhs"]  = rhs->name();
    map_fields_[rhs->name()] = rhs;
    m_assert(grid->IsAField(sol), "the solution field MUST exist in the grid");
    map_fields_[sol->name()] = sol;
    fields_nickname_["sol"]  = sol->name();
    m_assert(grid->IsAField(res), "the residual field MUST exist in the grid");
    map_fields_[res->name()] = res;
    fields_nickname_["res"]  = res->name();
    m_assert(rhs->lda() == sol->lda(), "the dimension between the source and the solution MUST match");
    m_assert(rhs->lda() == res->lda(), "the dimension between the source and the residual MUST match");

    // get by how many level I have to do
    p8est_t* forest      = grid->forest();
    sid_t    l_max_level = 0;
    for (p4est_topidx_t it = forest->first_local_tree; it <= forest->last_local_tree; it++) {
        // get the current tree
        p8est_tree_t* ctree = p8est_tree_array_index(forest->trees, it + forest->first_local_tree);
        // get the max delta level over the current tree
        l_max_level = m_max(ctree->maxlevel, l_max_level);  // max level is given by the tree
#ifndef NDEBUG
        // check that no quadrant is beneath the fft_level
        for (lid_t qid = 0; qid < ctree->quadrants.elem_count; qid++) {
            qdrt_t* quad = p8est_quadrant_array_index(&ctree->quadrants, qid);
            m_assert(fft_level_ <= quad->level, "the quad %d is beneath the level of FFT, not implemented yet", qid);
        }
#endif
    }
    // upate the counters on every rank
    m_assert(sizeof(l_max_level) == 1, "change the MPI datatype bellow");
    MPI_Allreduce(&l_max_level, &max_level_, 1, MPI_CHAR, MPI_MAX, MPI_COMM_WORLD);
    n_level_ = max_level_ - fft_level_;

    // allocate the grids and the partitionner
    grids_    = reinterpret_cast<Grid**>(m_calloc(sizeof(Grid*) * (n_level_ + 1)));
    families_ = reinterpret_cast<MGFamily**>(m_calloc(sizeof(MGFamily*) * n_level_));
    parts_    = reinterpret_cast<Partitioner**>(m_calloc(sizeof(Partitioner*) * n_level_));    

    // remember the original grid
    grids_[n_level_] = grid;

    // coarsen the grid
    for (sid_t il = (n_level_ - 1); il >= 0; il--) {
        // create a new grid and get the old information
        grids_[il] = new Grid();
        grids_[il]->CopyFrom(grids_[il + 1]);

        // regester myself in the user_pointer
        tmp_ilevel_               = il;
        p8est_t* curr_forest      = grids_[il]->forest();
        curr_forest->user_pointer = this;

        // get the number of CHILDREN on the fine level, i.e. the number of children entering the family
        sc_array_t quadarray = grids_[il+1]->mesh()->quad_level[fft_level + il + 1];
        m_assert(quadarray.elem_count < numeric_limits<lid_t>::max(), "the number of quad is too big");
        // create the new family
        families_[il] = new MGFamily(quadarray.elem_count);

        // coarsen by one level the quads on the highest level and only allocate the fields in the MG map
        p8est_coarsen_ext(curr_forest, 0, 0, cback_Level, nullptr, cback_MGCreateFamilly);
        m_assert(families_[il]->parent_count() == (quadarray.elem_count / P8EST_CHILDREN), "those two numbers must match: %d vs %ld", families_[il]->parent_count(), (quadarray.elem_count / P8EST_CHILDREN));

        // partition the grid and remember it, only the fields in the multigrid map exist now!
        parts_[il] = new Partitioner(&map_fields_, grids_[il], false);
        // create a consistent Ghost layout
        grids_[il]->SetupGhost();
        // the grid is now partitioned on a coarser level with new ghosts
    }
    //-------------------------------------------------------------------------
    m_end;
}

Multigrid::~Multigrid() {
    m_begin;
    //-------------------------------------------------------------------------

    for (sid_t il = (n_level_ - 1); il >= 0; il--) {
        delete (parts_[il]);
        delete (families_[il]);
        delete (grids_[il]);
    }

    m_free(parts_);
    m_free(families_);
    m_free(grids_);
    //-------------------------------------------------------------------------
    m_end;
}

void Multigrid::Solve() {
    m_begin;
    //-------------------------------------------------------------------------
    // get the fields
    Field* sol = map_fields_.at(fields_nickname_.at("sol"));
    Field* res = map_fields_.at(fields_nickname_.at("res"));
    Field* rhs = map_fields_.at(fields_nickname_.at("rhs"));
    // get which field will be send from one level to another
    map<string, Field*> map_sol_field;
    map_sol_field[fields_nickname_.at("sol")] = sol;
    map<string, Field*> map_rhs_field;
    map_rhs_field[fields_nickname_.at("rhs")] = rhs;

    Daxpy             daxpy_minus = Daxpy(-1.0);
    Daxpy             daxpy_plus  = Daxpy(+1.0);
    LaplacianCross<3> lapla       = LaplacianCross<3>();
    GaussSeidel<3>    gs          = GaussSeidel<3>(alpha_);

    IOH5 dump = IOH5("data");

    //---------------------------------------------------------------------------------------
    // downward pass
    for (sid_t il = n_level_; il > 0; il--) {
        Grid* grid = grids_[il];

        // IOH5 dump = IOH5("data");
        // string fname = "gs_" + std::to_string(il);
        // dump(grid,rhs,fname);

        ErrorCalculator error = ErrorCalculator();
        real_t          norm2;
        lapla(sol, res, grid);
        error.Norm2(grid, res, rhs, &norm2);
        m_log("entering level %d: error is %e", il, norm2);

        // Gauss-Seidel - laplacian(sol) = rhs
        for (sid_t ie = 0; ie < eta_1_; ie++) {
            gs(sol, rhs, nullptr, grid);
            // compute the error:
            lapla(sol, res, grid);  // laplacian(sol) = res
            error.Norm2(grid, res, rhs, &norm2);
            m_log("----> %e", norm2);
        }
        // compute the residual as rhs - A x
        lapla(sol, res, grid);             // laplacian(sol) = res
        dump(grid,res,"level1_down_lapla");
        daxpy_minus(grid, res, rhs, res);  // res = - 1.0 * res + rhs
        dump(grid,rhs,"level1_down_rhs");
        dump(grid,res,"level1_down_error");
        dump(grid,sol,"level1_down_sol");
        

        // do the interpolation
        families_[il - 1]->ToParents(res, rhs, grid->interp());

        // do the partitioning, send the new rhs
        parts_[il - 1]->Start(&map_rhs_field, M_FORWARD);
        parts_[il - 1]->End(&map_rhs_field, M_FORWARD);
    }
    //---------------------------------------------------------------------------------------
    // do the direct solve
    Grid*           grid_fft = grids_[0];
    ErrorCalculator error2   = ErrorCalculator();
    real_t          norm22;
    lapla(sol, res, grid_fft);
    error2.Norm2(grid_fft, res, rhs, &norm22);
    m_log("entering level %d: error is %e", 0, norm22);

    dump(grid_fft,sol,"fft_before_sol");
    dump(grid_fft,res,"fft_before_lapla");
    dump(grid_fft,rhs,"fft_before_rhs");
    // Gauss-Seidel - laplacian(sol) = rhs
    for (sid_t ie = 0; ie < 5; ie++) {
        ErrorCalculator error = ErrorCalculator();
        real_t          norm2;
        // lapla(sol, res, grids_[0]);
        // error.Norm2(grid_fft, res, rhs, &norm2);
        // m_log("error at level %d before is %e", 0, norm2);
        gs(sol, rhs,nullptr, grid_fft);
        // IOH5   dump  = IOH5("data");
        // string fname = "sol_gs_" + std::to_string(il) + "_" + std::to_string(ie);
        // dump(grid, sol, fname);
        // compute the error:
        lapla(sol, res, grid_fft);  // laplacian(sol) = res
        // ErrorCalculator error = ErrorCalculator();
        // real_t          norm2;
        error.Norm2(grid_fft, res, rhs, &norm2);
        
        m_log("------> %e", norm2);
    }
    dump(grid_fft,sol,"fft_sol");
    dump(grid_fft,res,"fft_lapla");
    dump(grid_fft,rhs,"fft_rhs");

    //---------------------------------------------------------------------------------------
    // downward pass
    for (sid_t il = 1; il <= n_level_; il++) {
        // get the ghost for the solution as we need them to interpolate!!
        grids_[il-1]->GhostPull(sol);
         // do the partitioning, send the solution + ghost
        parts_[il - 1]->Start(&map_sol_field, M_BACKWARD);
        parts_[il - 1]->End(&map_sol_field, M_BACKWARD);

        // do the interpolation
        Grid* grid = grids_[il];
        families_[il - 1]->ToChildren(sol, res, grid->interp());

        dump(grid,res,"solinterp");

        ErrorCalculator error = ErrorCalculator();
        real_t          norm2;
        // lapla(sol, res, grid);
        // error.Norm2(grid, res, rhs, &norm2);
        // m_log("entering level %d: error is %e", il, norm2);

        // update the solution with it's coarsegrid version
        daxpy_plus(grid, res, sol, sol);  // sol = res + sol

        dump(grid,sol,"level1_up_sol");
        dump(grid,res,"level1_up_solinterp");

        lapla(sol, res, grid);
        error.Norm2(grid, res, rhs, &norm2);
        dump(grid, res, "level1_up_error_lapla");
        dump(grid, rhs, "level1_up_error_rhs");
        m_log("after correction error is %e", norm2);

        // Gauss-Seidel - laplacian(sol) = rhs
        for (sid_t ie = 0; ie < eta_2_; ie++) {
            gs(sol, rhs,nullptr, grid);
            // compute the error:
            lapla(sol, res, grid);  // laplacian(sol) = res
            error.Norm2(grid, res, rhs, &norm2);
            // m_log("error at level %d after is %e", il, norm2);
            m_log("------> %e", norm2);
        }
    }
    //-------------------------------------------------------------------------
    m_end;
}