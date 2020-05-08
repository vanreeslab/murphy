#include "multigrid.hpp"

#include <mpi.h>

#include "defs.hpp"
#include "gridcallback.hpp"
#include "laplacian.hpp"
#include "daxpy.hpp"
#include "gaussseidel.hpp"

Multigrid::Multigrid(Grid* grid, sid_t fft_level,Field* src, Field* sol,Field* res) {
    m_begin;
    //-------------------------------------------------------------------------
    // store the desired fft level and which field will be used as what
    fft_level_        = fft_level;
    m_assert(grid->IsAField(src),"the source field MUST exist in the grid");
    map_fields_["src"] = src;
     m_assert(grid->IsAField(sol),"the solution field MUST exist in the grid");
    map_fields_["sol"] = sol;
    m_assert(grid->IsAField(res),"the residual field MUST exist in the grid");
    map_fields_["res"] = res;
    m_assert(src->lda() == sol->lda(),"the dimension between the source and the solution MUST match");
    m_assert(src->lda() == res->lda(),"the dimension between the source and the residual MUST match");
    
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
        // create a new grid with the link to the previous blocks
        grids_[il] = new Grid(grids_[il + 1]);
        // regester myself in the user_pointer
        tmp_ilevel_               = il;
        p8est_t* curr_forest      = grids_[il]->forest();
        curr_forest->user_pointer = this;
        // create a new family knowing the number of CHILDREN on the fine level
        sc_array_t quadarray = grid->mesh()->quad_level[fft_level + il + 1];
        families_[il] = new MGFamily(quadarray.elem_count);
        // coarsen by one level the quads on the highest level and only allocate the fields in the MG map
        p8est_coarsen_ext(curr_forest,0,0,cback_Level,nullptr,cback_MGCreateFamilly);
        m_assert(families_[il]->parent_count() == quadarray.elem_count,"those two numbers must match");
        // partition the grid and remember it, only the fields in the multigrid map exist now!
        parts_[il] = new Partitioner(&map_fields_,grids_[il],false);
        // create a consistent Ghost layout


        // the grid is now partitioned on a coarser level!
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
    Field* sol = map_fields_.at("sol");
    Field* res = map_fields_.at("res");
    Field* src = map_fields_.at("src");
    // get which field will be send from one level to another
    map<string,Field*> part_fields;
    part_fields["res"] = res;

    Daxpy             daxpy = Daxpy(-1.0);
    LaplacianCross<5> lapla = LaplacianCross<5>();
    GaussSeidel<5> gs = GaussSeidel<5>(1.95);

    // get the initial residual

    lapla(sol,res,grids_[n_level_]);
    daxpy(grids_[n_level_],sol,src,res);

    // downward pass
    for (sid_t il = (n_level_ - 1); il >= 0; il--) {
        Grid* grid = grids_[il];
        
        // compute some GS on the solution
        for(sid_t ie = 0; ie< eta_1_; ie++){
            gs(res,sol,grid);
        }

        // compute the residual as b - A x
        lapla(sol,res,grid);
        daxpy(grid,sol,src,res);

        // do the interpolation
        families_[il]->ToParents(res, grid->interp());
        // do the partitioning
        parts_[il]->Start(&part_fields);
        parts_[il]->End(&part_fields);
    }
    // do the direct solve

    // upward pass
    for (sid_t il = 0; il < n_level_; il--) {
        Grid* grid = grids_[il];
        
        // compute some GS on the solution
        for(sid_t ie = 0; ie< eta_1_; ie++){
            gs(res,sol,grid);
        }

        // compute the residual as b - A x
        lapla(sol,res,grid);
        daxpy(grid,sol,src,res);

        // do the interpolation
        families_[il]->ToParents(res, grid->interp());
        // do the partitioning
        parts_[il]->Start(&part_fields);
        parts_[il]->End(&part_fields);
    }
    //-------------------------------------------------------------------------
    m_end;
}