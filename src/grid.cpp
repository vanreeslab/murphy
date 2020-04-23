#include "grid.hpp"

#include <omp.h>
#include <p8est_extended.h>

#include "gridcallback.hpp"
#include "operator.hpp"
#include "wavelet.hpp"
#include "partitioner.hpp"

using std::string;
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
Grid::Grid(const lid_t ilvl, const bool isper[3], const lid_t l[3], MPI_Comm comm, Prof* prof)
    : ForestGrid(ilvl, isper, l, sizeof(GridBlock*), comm) {
    m_begin;
    //-------------------------------------------------------------------------
    // create a default interpolator
    interp_ = new Wavelet<3>();

    // profiler
    prof_ = prof;

    // create the associated blocks
    p8est_iterate(forest_, NULL, NULL, cback_CreateBlock, NULL, NULL, NULL);

    // create the ghosts structure
    ghost_ = new Ghost(this);
    //-------------------------------------------------------------------------
    m_log("uniform grid created with %ld blocks on %ld trees using %d ranks and %d threads", forest_->global_num_quadrants, forest_->trees->elem_count, forest_->mpisize, omp_get_max_threads());
    m_end;
}

/**
 * @brief Destroy the grid_t object: destroy the p8est objects
 * 
 */
Grid::~Grid() {
    m_begin;
    //-------------------------------------------------------------------------
    // destroy the interpolator
    delete (interp_);
    delete (ghost_);
    // destroy the remaining blocks
    p8est_iterate(forest_, NULL, NULL, cback_DestroyBlock, NULL, NULL, NULL);
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief returns the size of the memory (in bytes) taken by the grid
 * 
 * @return size_t 
 */
size_t Grid::LocalMemSize() const {
    m_begin;
    //-------------------------------------------------------------------------
    size_t memsize = 0;

    memsize += sizeof(fields_);
    memsize += sizeof(prof_);
    // memsize += ghost_->LocalMemSize();
    // memsize += interp_->LocalMemSize();
    for (auto fid = fields_.begin(); fid != fields_.end(); fid++) {
        memsize += forest_->local_num_quadrants * (M_N * M_N * M_N) * fid->second->lda() * sizeof(real_t);
    }
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
    return (fields_.find(key) != fields_.end());
    // return (std::find(fields_.begin(),fields_.end(),field) != fields_.end());
}

void Grid::AddField(Field* field) {
    m_begin;
    //-------------------------------------------------------------------------
    if (!IsAField(field)) {
        // add the field
        fields_[field->name()] = field;
        // allocate the field on every block
        LoopOnGridBlock_(&GridBlock::AddField, field);
        m_verb("field %s has been added to the grid", field->name().c_str());
    } else {
        m_verb("field %s is already in the grid", field->name().c_str());
    }
    //-------------------------------------------------------------------------
    m_end;
}

void Grid::DeleteField(Field* field) {
    m_begin;
    //-------------------------------------------------------------------------
    if (IsAField(field)) {
        // add the field to everyblock
        LoopOnGridBlock_(&GridBlock::DeleteField, field);
        // remove the field
        fields_.erase(field->name());
    } else {
        m_verb("field %s is not in the grid", field->name().c_str());
    }
    //-------------------------------------------------------------------------
    m_end;
}

void Grid::GhostPull(Field* field) {
    m_begin;
    m_assert(interp_ != nullptr, "the inteprolator cannot be null");
    m_assert(IsAField(field), "the field does not belong to this grid");
    //-------------------------------------------------------------------------
    // if already computed, return
    if (field->ghost_status() || field->lda() < 1) {
        m_log("field %s has already valid ghosts", field->name().c_str());
        return;
    } else {
        // m_log("pulling ghosts for field %s", field->name().c_str());
        // start the send of the first dimension
        ghost_->PushToMirror(field, 0);
        ghost_->MirrorToGhostSend();

        for (int ida = 1; ida < field->lda(); ida++) {
            // receive the current communication, the mirrors are now free
            ghost_->MirrorToGhostRecv();
            // fill the mirror and initiate the next send
            ghost_->PushToMirror(field, ida);
            ghost_->MirrorToGhostSend();
            // handle the last dimension just received
            ghost_->PullFromGhost(field, ida - 1, interp_);
        }
        // receive and end the last dimension
        ghost_->MirrorToGhostRecv();
        ghost_->PullFromGhost(field, field->lda() - 1, interp_);

        // set that everything is ready for the field
        field->ghost_status(true);
        // m_log("ghosts pulling is done for field %s", field->name().c_str());
    }
    //-------------------------------------------------------------------------
    m_end;
}

void Grid::Refine(const sid_t delta_level) {
    m_begin;
    //-------------------------------------------------------------------------
    // we create the new blocks
    for (int id = 0; id < delta_level; id++) {
        // compute the ghost needed by the interpolation
        for (auto fid = fields_.begin(); fid != fields_.end(); fid++) {
            // set the working field before entering the callback
            // working_callback_field_ = fid->second;
            // get the ghosts needed by the interpolation
            GhostPull(fid->second);
        }
        // delete the soon-to be outdated ghost and mesh
        delete (ghost_);
        ResetP4estGhostMesh();
        // set the grid in the forest for the callback
        forest_->user_pointer = (void*)this;
        // do the p4est interpolation by callback
        p8est_refine_ext(forest_, 0, P8EST_MAXLEVEL, cback_Yes, nullptr, cback_Interpolate);
        // balance the partition
        p8est_balance_ext(forest_, P8EST_CONNECT_FULL, NULL, cback_Interpolate);
        // partition the grid
        Partitioner partition(&fields_, this);
        partition.Start(&fields_);
        partition.End(&fields_);
        // create a new ghost and mesh
        SetupP4estGhostMesh();
        ghost_ = new Ghost(this);
        // set the ghosting fields as non-valid
        for (auto fid = fields_.begin(); fid != fields_.end(); fid++) {
            fid->second->ghost_status(false);
        }
    }
    //-------------------------------------------------------------------------
    m_log("refined grid created with %ld blocks on %ld trees using %d ranks and %d threads", forest_->global_num_quadrants, forest_->trees->elem_count, forest_->mpisize, omp_get_max_threads());
    m_end;
}

void Grid::Coarsen(const sid_t delta_level) {
    m_begin;
    //-------------------------------------------------------------------------
    // we create the new blocks
    for (int id = 0; id < delta_level; id++) {
        // compute the ghost needed by the interpolation
        for (auto fid = fields_.begin(); fid != fields_.end(); fid++) {
            // set the working field before entering the callback
            // working_callback_field_ = fid->second;
            // get the ghosts needed by the interpolation
            GhostPull(fid->second);
        }
        // delete the soon-to be outdated ghost and mesh
        delete (ghost_);
        ResetP4estGhostMesh();
        // set the grid in the forest for the callback
        forest_->user_pointer = (void*)this;
        // do the p4est interpolation by callback
        p8est_coarsen_ext(forest_, 0, 0, cback_Yes, nullptr, cback_Interpolate);
        // balance the partition
        p8est_balance_ext(forest_, P8EST_CONNECT_FULL, NULL, cback_Interpolate);
        // partition the grid
        Partitioner partition(&fields_, this);
        partition.Start(&fields_);
        partition.End(&fields_);
        // create a new ghost and mesh
        SetupP4estGhostMesh();
        ghost_ = new Ghost(this);
        // set the ghosting fields as non-valid
        for (auto fid = fields_.begin(); fid != fields_.end(); fid++) {
            fid->second->ghost_status(false);
        }
    }
    //-------------------------------------------------------------------------
    m_log("coarsened grid created with %ld blocks on %ld trees using %d ranks and %d threads", forest_->global_num_quadrants, forest_->trees->elem_count, forest_->mpisize, omp_get_max_threads());
    m_end;
}

void Grid::SetTol(const real_t refine_tol, const real_t coarsen_tol) {
    m_begin;
    //-------------------------------------------------------------------------
    rtol_ = refine_tol;
    ctol_ = coarsen_tol;
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief adapt = (refine or coarsen once) each block, based on the given field and the tolerance
 * 
 * @param field 
 * @param tol 
 */
void Grid::Adapt(Field* field) {
    m_begin;
    m_log("grid adaptation started...");
    //-------------------------------------------------------------------------
    // store the criterion field
    tmp_field_ = field;
    // compute the ghost needed by the interpolation of everyblock
    for (auto fid = fields_.begin(); fid != fields_.end(); fid++) {
        GhostPull(fid->second);
    }
    // delete the soon-to be outdated ghost and mesh
    delete (ghost_);
    ResetP4estGhostMesh();
    // set the grid in the forest for the callback
    forest_->user_pointer = (void*)this;
    // coarsen the needed block
    p8est_coarsen_ext(forest_, 0, 0, cback_Wavelet, nullptr, cback_Interpolate);
    // refine the needed blocks
    p8est_refine_ext(forest_, 0, P8EST_MAXLEVEL, cback_Wavelet, nullptr, cback_Interpolate);
    // balance the partition
    p8est_balance_ext(forest_, P8EST_CONNECT_FULL, NULL, cback_Interpolate);
    // partition the grid
    Partitioner partition(&fields_, this);
    partition.Start(&fields_);
    partition.End(&fields_);
    // create a new ghost and mesh
    SetupP4estGhostMesh();
    ghost_ = new Ghost(this);
    // set the ghosting fields as non-valid
    for (auto fid = fields_.begin(); fid != fields_.end(); fid++) {
        fid->second->ghost_status(false);
    }
    //-------------------------------------------------------------------------
    m_log("...grid adaptation done: now %ld blocks on %ld trees using %d ranks and %d threads", forest_->global_num_quadrants, forest_->trees->elem_count, forest_->mpisize, omp_get_max_threads());
    m_end;
}

/**
 * @brief iterates on the blocks and performs a simple @ref bop_t operation
 * 
 * @warning for allocation and block management only. Use Operators for computations
 * 
 * @param op 
 * @param grid 
 * @param field 
 */
void Grid::LoopOnGridBlock_(const bop_t op, Field* field) const {
    m_begin;
    //-------------------------------------------------------------------------
    // get the grid info
    for (p4est_topidx_t it = forest_->first_local_tree; it <= forest_->last_local_tree; it++) {
        p8est_tree_t* tree    = p8est_tree_array_index(forest_->trees, it);
        const size_t  nqlocal = tree->quadrants.elem_count;

#pragma omp parallel for firstprivate(field)
        for (size_t bid = 0; bid < nqlocal; bid++) {
            p8est_quadrant_t* quad  = p8est_quadrant_array_index(&tree->quadrants, bid);
            GridBlock*        block = *(reinterpret_cast<GridBlock**>(quad->p.user_data));
            // apply
            (block->*op)(field);
        }
        // downgrade the ghost status since we changed its value
        field->ghost_status(false);
        //-------------------------------------------------------------------------
        m_end;
    }
}