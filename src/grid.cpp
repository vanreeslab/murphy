#include "grid.hpp"

#include <p8est_extended.h>
#include <omp.h>

#include "block.hpp"
#include "operator.hpp"

using std::string;

void cback_CreateBlock(p8est_iter_volume_info_t * info,void *user_data) {
    m_begin;
    //-------------------------------------------------------------------------
    p8est_t*              forest     = info->p4est;
    p8est_quadrant_t*     quad       = info->quad;
    p4est_topidx_t        which_tree = info->treeid;
    p8est_connectivity_t* connect    = forest->connectivity;
    // get the starting position
    real_t xyz[3];
    p8est_qcoord_to_vertex(connect, which_tree, quad->x, quad->y, quad->z, xyz);

    real_t len        = P8EST_QUADRANT_LEN(quad->level)*(1.0/P8EST_ROOT_LEN);
    quad->p.user_data = new Block(len, xyz, quad->level);
    //-------------------------------------------------------------------------
    m_end;
}

void cback_DestroyBlock(p8est_iter_volume_info_t* info, void* user_data) {
    m_begin;
    //-------------------------------------------------------------------------
    p8est_quadrant_t* quad = info->quad;
    delete ((Block*)quad->p.user_data);
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief Construct a new grid_t: initialize the p8est objects
 * 
 * The grid is initialized at @ref ilvl as a uniform grid.
 * It means that the number of blocks = (l[0]*l[1]*l[2]) * (2^ilvl)^3
 * 
 * @param ilvl the initialization level, for every tree
 * @param isper isper[i] indicates that the ith direction is periodic (x:0 y:1 z:2)
 * @param L the number of trees in each direction, i.e. the aspect ratio of the domain
 * @param comm the communicator to use
 * @param prof an existing profiler for the programm
 */
Grid::Grid(const lid_t ilvl, const bool isper[3], const lid_t l[3], MPI_Comm comm, Prof* prof) {
    m_begin;
    m_assert(ilvl >= 0, "the init level has to be >= 0");
    m_assert(ilvl <= P8EST_MAXLEVEL, "the init level has to be <= P8EST_MAXLEVEL(=19)");
    //-------------------------------------------------------------------------
    prof_ = prof;

    // create the connect as a box of L[0]xL[1]xL[2] trees
    p8est_connectivity_t* connect = p8est_connectivity_new_brick(l[0], l[1], l[2], isper[0], isper[1], isper[2]);

    // create the forest at a given level, the associated ghost and mesh object
    forest_ = p8est_new_ext(comm, connect, 0, ilvl, 1, sizeof(Block*), nullptr, nullptr);
    ghost_  = p8est_ghost_new(forest_, P8EST_CONNECT_FULL);
    mesh_   = p8est_mesh_new_ext(forest_, ghost_, 1, 1, P8EST_CONNECT_FULL);

    // declare the mesh as valid
    is_mesh_valid_ = true;

    // create the associated blocks
    p8est_iterate(forest_, NULL, NULL, cback_CreateBlock, NULL, NULL, NULL);

    //-------------------------------------------------------------------------
    m_log("uniform grid created with %ld blocks on %ld trees using %d ranks and %d threads", forest_->global_num_quadrants, forest_->trees->elem_count, forest_->mpisize,omp_get_max_threads());
    m_end;
}

/**
 * @brief Destroy the grid_t object: destroy the p8est objects
 * 
 */
Grid::~Grid() {
    m_begin;
    //-------------------------------------------------------------------------
    // destroy the remaining blocks
    p8est_iterate(forest_, NULL, NULL, cback_DestroyBlock, NULL, NULL, NULL);
    // destroy the structures
    p8est_mesh_destroy(mesh_);
    p8est_ghost_destroy(ghost_);
    // destroy the connectivity and the forest
    p8est_connectivity_t* connect = forest_->connectivity;
    p8est_destroy(forest_);
    p8est_connectivity_destroy(connect);
    //-------------------------------------------------------------------------
    m_end;
}

size_t Grid::LocalMemSize() const {
    m_begin;
    //-------------------------------------------------------------------------
    size_t memsize = 0;
    memsize += sizeof(forest_);
    memsize += sizeof(ghost_);
    memsize += sizeof(mesh_);
    memsize += sizeof(prof_);
    memsize += sizeof(lda_) + lda_.size() * sizeof(sid_t);
    memsize += forest_->local_num_quadrants * (M_N * M_N * M_N) * sizeof(real_t);
    //-------------------------------------------------------------------------
    m_end;
    return memsize;
}

size_t Grid::LocalNumDof() const {
    return forest_->local_num_quadrants * (M_N * M_N * M_N);
}
size_t Grid::GlobalNumDof() const {
    return forest_->global_num_quadrants * (M_N * M_N * M_N);
}

bool Grid::IsAField(const Field* field) const {
    std::string key = field->name();
    return (lda_.find(key) != lda_.end());
}

void Grid::AddField(Field* field) {
    m_begin;
    //-------------------------------------------------------------------------
    if (!IsAField(field)) {
        //get the key = name of the field
        string key = field->name();
        // create a new lda entry
        lda_[key] = field->lda();
        // add the field to everyblock
        DoOp<nullptr_t>(&Block::AddField, this, field, nullptr);
    }
    //-------------------------------------------------------------------------
    m_end;
}

void Grid::DeleteField(Field* field) {
    m_begin;
    //-------------------------------------------------------------------------
    if (IsAField(field)) {
        //get the key = name of the field
        string key = field->name();
        // create a new lda entry
        lda_.erase(key);
        // add the field to everyblock
        DoOp<nullptr_t>(&Block::DeleteField, this, field, nullptr);
    }
    //-------------------------------------------------------------------------
    m_end;
}