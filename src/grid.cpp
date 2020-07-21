#include "grid.hpp"

#include <omp.h>
#include <p8est_extended.h>

#include "gridcallback.hpp"
#include "partitioner.hpp"
#include "wavelet.hpp"

using std::string;

/**
 * @brief Construct a new Grid which is empty, only the interpolators have been associated
 * 
 */
Grid::Grid() : ForestGrid() {
    //-------------------------------------------------------------------------
    prof_     = nullptr;
    ghost_    = nullptr;
    // create a default interpolator
    interp_ = new Wavelet<2, 2>();
    //-------------------------------------------------------------------------
};

/**
 * @brief Construct a new Grid a a uniform grid, distributed among the cpus
 * 
 * Initialize the ForestGrid, the interpolator, the profiler, the ghosts and the grid blocks.
 * The grid is defined as a set of l[0]xl[1]xl[2] octrees, each of them refined up to level ilvl.
 * 
 * @param ilvl the starting initialization level
 * @param isper defines the periodicity of the dimension (cannot be changed afterwards)
 * @param l the aspect ratio of the domain (integers), defining how many trees are used in each direction
 * @param comm the MPI communicator used
 * @param prof the profiler pointer if any (can be nullptr)
 */
Grid::Grid(const lid_t ilvl, const bool isper[3], const lid_t l[3], MPI_Comm comm, Prof* prof)
    : ForestGrid(ilvl, isper, l, sizeof(GridBlock*), comm) {
    m_begin;
    //-------------------------------------------------------------------------
    // profiler
    prof_ = prof;
    // init the profiler tracking
    if (prof_ != nullptr) {
        // ghost
        prof->Create("ghost_init");
        prof->Create("ghost_comm");
        prof->Create("ghost_cmpt");
        prof->Create("ghost_comm_start", "ghost_comm");
        prof->Create("ghost_comm_wait", "ghost_comm");
        // p4est calls
        prof->Create("p4est_refcoarse");
        prof->Create("p4est_balance");
        prof->Create("p4est_partition_init");
        prof->Create("p4est_partition_comm");
        prof->Create("cback_interpolate","p4est_refcoarse");
        // stencils
        prof->Create("stencil_inner");
        prof->Create("stencil_outer");
    }
    // create a default interpolator
    interp_ = new Wavelet<2,2>();
    // create the associated blocks
    p8est_iterate(forest_, NULL, NULL, cback_CreateBlock, NULL, NULL, NULL);
    // partition the grid to have compatible grid
    Partitioner part = Partitioner(&fields_,this,true);
    part.Start(&fields_,M_FORWARD);
    part.End(&fields_,M_FORWARD);
    // setup the ghost stuctures as the mesh will not change anymore
    SetupGhost();
    //-------------------------------------------------------------------------
    m_log("uniform grid created with %ld blocks on %ld trees using %d ranks and %d threads", forest_->global_num_quadrants, forest_->trees->elem_count, forest_->mpisize, omp_get_max_threads());
    m_end;
}

/**
 * @brief Copy the ForestGrid part from a grid
 * 
 * @param grid the source grid
 */
void Grid::CopyFrom(Grid* grid){
    m_begin;
    //-------------------------------------------------------------------------
    this->ForestGrid::CopyFrom(grid);
    // copy the field mapping
    for (auto iter = grid->FieldBegin(); iter != grid->FieldEnd(); iter++) {
        string name   = iter->first;
        Field* fid    = iter->second;
        fields_[name] = fid;
    }
    // copy the profiler
    prof_ = grid->profiler();
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief Destroy the Grid, frees all the blocks and the fields contained (if not done)
 * 
 * @warning if we are a copy of another grid, we cannot know which block has been allocated by us,
 * so we do not free the blocks!!
 * 
 */
Grid::~Grid() {
    m_begin;
    //-------------------------------------------------------------------------
    // destroy the interpolator and the details they are mine
    if (interp_ != nullptr) {
        m_log("dealloc the interp");
        delete (interp_);
    }
    // destroy the ghosts, they are mine as well
    DestroyGhost();
    // destroy the remaining blocks
    if (is_connect_owned_) {
        p8est_iterate(forest_, NULL, NULL, cback_DestroyBlock, NULL, NULL, NULL);
    }
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief setup the Ghost structure when the mesh is not going to change anymore
 * 
 * @warning this function cannot be called on on existing structure
 * 
 */
void Grid::SetupGhost() {
    m_begin;
    m_assert(ghost_ == nullptr, "cannot create something that already exists");
    //-------------------------------------------------------------------------
    // create the forestGrid part
    this->SetupP4estGhostMesh();
    // create the ghosts structure
    m_profStart(prof_,"ghost_init");
    ghost_ = new Ghost(this, interp_);
    m_profStop(prof_,"ghost_init");
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief destroys the Ghost structure of both the ForestGrid and the grid
 * 
 */
void Grid::DestroyGhost() {
    m_begin;
    //-------------------------------------------------------------------------
    if (ghost_ != nullptr) {
        m_log("dealloc the ghost");
        delete (ghost_);
        ghost_ = nullptr;
    }
    this->ResetP4estGhostMesh();
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief returns the local memory size (in Byte) of the grid, including every field contained in it
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

/**
 * @brief returns the number of local number of degree of freedom for one field containing one dimension
 * 
 * @return size_t 
 */
size_t Grid::LocalNumDof() const {
    return forest_->local_num_quadrants * (M_N * M_N * M_N);
}

/**
 * @brief returns the number of global number of degree of freedom for one field containing one dimension
 * 
 * @return size_t 
 */
size_t Grid::GlobalNumDof() const {
    return forest_->global_num_quadrants * (M_N * M_N * M_N);
}

bool Grid::IsAField(const Field* field) const {
    std::string key = field->name();
    return (fields_.find(key) != fields_.end());
    // return (std::find(fields_.begin(),fields_.end(),field) != fields_.end());
}

/**
 * @brief register a field in the current grid
 *
 * It loops on the blocks to allocate the needed memory
 * 
 * @param field 
 */
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

/**
 * @brief unregister a field in the current grid and free the associated memory on each block
 * 
 * @param field
 */
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

/**
 * @brief reset the fields_ list to the given list of new fields
 * 
 * @warning if mis-used, this function will break the whole run...
 * 
 * @param fields 
 */
void Grid::ResetFields(const map<string, Field*>* fields) {
    m_begin;
    //-------------------------------------------------------------------------
    // clear the current map
    fields_.clear();
    // copy the new one
    for (auto iter = fields->begin(); iter != fields->end(); iter++) {
        fields_[iter->first] = iter->second;
        // check if we satisfy the requirements on the key
        m_assert(iter->first == iter->second->name(),"the key of the map must be the name of the field");
    }
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief Pull ghost points (take the values from the neighbors): start to obtain the ghost points from coarser and same level neighbors
 * 
 * @warning you have to call the @ref GhostPull_CoarseEnd() before making any other communication
 * 
 * @param field the considered field
 * @param ida the dimension of the field which has to be send
 */
void Grid::GhostPull_CoarseStart(Field* field, const sid_t ida) {
    m_begin;
    m_assert(0 <= ida && ida < field->lda(), "the ida is not within the field's limit");
    m_assert(field != nullptr, "the source field cannot be null");
    m_assert(IsAField(field), "the field does not belong to this grid");
    m_assert(ghost_ != nullptr, "The ghost structure is not valid, unable to use it");
    //-------------------------------------------------------------------------
    if (!field->ghost_status()) {
        m_profStart(prof_,"ghostpull_coarse_start");
        ghost_->GetGhost_Post(field,ida);
        m_profStop(prof_,"ghostpull_coarse_start");
    }
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief Pull ghost points (take the values from the neighbors): end to obtain the ghost points from coarser and same level neighbors
 * 
 * @param field the considered field
 * @param ida the received dimension
 */
void Grid::GhostPull_CoarseEnd(Field* field, const sid_t ida) {
    m_begin;
    m_assert(0 <= ida && ida < field->lda(), "the ida is not within the field's limit");
    m_assert(field != nullptr, "the source field cannot be null");
    m_assert(IsAField(field), "the field does not belong to this grid");
    m_assert(ghost_ != nullptr, "The ghost structure is not valid, unable to use it");
    //-------------------------------------------------------------------------
    if (!field->ghost_status()) {
        m_profStart(prof_, "ghostpull_coarse_end");
        ghost_->GetGhost_Wait(field, ida);
        m_profStop(prof_, "ghostpull_coarse_end");
    }
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief Pull ghost points (take the values from the neighbors): start to obtain the ghost points from finer neighbors
 * 
 * @param field the considered field
 * @param ida the received dimension
 */
void Grid::GhostPull_FineStart(Field* field, const sid_t ida) {
    m_begin;
    m_assert(0 <= ida && ida < field->lda(), "the ida is not within the field's limit");
    m_assert(field != nullptr, "the source field cannot be null");
    m_assert(IsAField(field), "the field does not belong to this grid");
    m_assert(ghost_ != nullptr, "The ghost structure is not valid, unable to use it");
    //-------------------------------------------------------------------------
    if (!field->ghost_status()) {
        m_profStart(prof_, "ghostpull_fine_start");
        ghost_->PutGhost_Post(field, ida);
        m_profStop(prof_, "ghostpull_fine_start");
    }
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief Pull ghost points (take the values from the neighbors): end to obtain the ghost points from finer neighbors
 * 
 * @param field the considered field
 * @param ida the filled dimension
 */
void Grid::GhostPull_FineEnd(Field* field, const sid_t ida) {
    m_begin;
    m_assert(0 <= ida && ida < field->lda(), "the ida is not within the field's limit");
    m_assert(field != nullptr, "the source field cannot be null");
    m_assert(IsAField(field), "the field does not belong to this grid");
    m_assert(ghost_ != nullptr, "The ghost structure is not valid, unable to use it");
    //-------------------------------------------------------------------------
    if (!field->ghost_status()) {
        m_profStart(prof_, "ghostpull_fine_start");
        ghost_->PutGhost_Wait(field, ida);
        m_profStop(prof_, "ghostpull_fine_start");
    }
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief Pull the ghost points (take the values from the neighbors): after this function, the ghost status of the field is set as up-to-date
 * 
 * @param field the field which requires the ghost
 */
void Grid::GhostPull(Field* field) {
    m_begin;
    m_assert(field != nullptr, "the source field cannot be null");
    m_assert(ghost_ != nullptr, "The ghost structure is not valid, unable to use it");
    //-------------------------------------------------------------------------
    // start the send in the first dimension
    for (int ida = 0; ida < field->lda(); ida++) {
        GhostPull_CoarseStart(field, ida);
        GhostPull_CoarseEnd(field, ida);
        GhostPull_FineStart(field, ida);
        GhostPull_FineEnd(field, ida);
    }
    // // set that everything is ready for the field
    field->ghost_status(true);
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief Refine the grid by delta_level locally
 * 
 * @param delta_level the number of level every block will be refined
 */
void Grid::Refine(const sid_t delta_level) {
    m_begin;
    //-------------------------------------------------------------------------
    // we create the new blocks
    for (int id = 0; id < delta_level; id++) {
        // compute the ghost needed by the interpolation
        for (auto fid = fields_.begin(); fid != fields_.end(); fid++) {
            GhostPull(fid->second);
        }
        // delete the soon-to be outdated ghost and mesh
        DestroyGhost();
        // set the grid in the forest for the callback
        forest_->user_pointer = reinterpret_cast<void*>(this);
        // do the p4est interpolation by callback
        if (prof_ != nullptr) {
            prof_->Start("p4est_refcoarse");
        }
        p8est_refine_ext(forest_, 0, P8EST_MAXLEVEL, cback_Yes, nullptr, cback_Interpolate);
        if (prof_ != nullptr) {
            prof_->Stop("p4est_refcoarse");
            prof_->Start("p4est_balance");
        }
        // balance the partition
        p8est_balance_ext(forest_, P8EST_CONNECT_FULL, NULL, cback_Interpolate);
        if (prof_ != nullptr) {
            prof_->Stop("p4est_balance");
            prof_->Start("p4est_partition_init");
        }
        // partition the grid
        Partitioner partition(&fields_, this, true);
        if (prof_ != nullptr) {
            prof_->Stop("p4est_partition_init");
            prof_->Start("p4est_partition_comm");
        }
        partition.Start(&fields_,M_FORWARD);
        partition.End(&fields_,M_FORWARD);
        if (prof_ != nullptr) {
            prof_->Stop("p4est_partition_comm");
        }
        // create a new ghost and mesh
        SetupGhost();
        // set the ghosting fields as non-valid
        for (auto fid = fields_.begin(); fid != fields_.end(); fid++) {
            fid->second->ghost_status(false);
        }
    }
    //-------------------------------------------------------------------------
    m_log("refined grid created with %ld blocks on %ld trees using %d ranks and %d threads", forest_->global_num_quadrants, forest_->trees->elem_count, forest_->mpisize, omp_get_max_threads());
    m_end;
}

/**
 * @brief coarsen the grid by delta_level locally
 * 
 * @param delta_level 
 */
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
        DestroyGhost();
        // set the grid in the forest for the callback
        forest_->user_pointer = reinterpret_cast<void*>(this);
        // do the p4est interpolation by callback
        if (prof_ != nullptr) {
            prof_->Start("p4est_refcoarse");
        }
        p8est_coarsen_ext(forest_, 0, 0, cback_Yes, nullptr, cback_Interpolate);
        if (prof_ != nullptr) {
            prof_->Stop("p4est_refcoarse");
            prof_->Start("p4est_balance");
        }
        // balance the partition
        p8est_balance_ext(forest_, P8EST_CONNECT_FULL, NULL, cback_Interpolate);
        if (prof_ != nullptr) {
            prof_->Stop("p4est_balance");
            prof_->Start("p4est_partition_init");
        }
        // partition the grid
        Partitioner partition(&fields_, this, true);
        if (prof_ != nullptr) {
            prof_->Stop("p4est_partition_init");
            prof_->Start("p4est_partition_comm");
        }
        partition.Start(&fields_,M_FORWARD);
        partition.End(&fields_,M_FORWARD);
        // create a new ghost and mesh
        if (prof_ != nullptr) {
            prof_->Stop("p4est_partition_comm");
        }
        SetupGhost();
        // set the ghosting fields as non-valid
        for (auto fid = fields_.begin(); fid != fields_.end(); fid++) {
            fid->second->ghost_status(false);
        }
    }
    //-------------------------------------------------------------------------
    m_log("coarsened grid created with %ld blocks on %ld trees using %d ranks and %d threads", forest_->global_num_quadrants, forest_->trees->elem_count, forest_->mpisize, omp_get_max_threads());
    m_end;
}

/**
 * @brief sets the refinement tolerance for grid adaptation (see @ref Adapt)
 * 
 * @param refine_tol we refine if criterion > refine_tol
 * @param coarsen_tol we coarsen if criterion < coarsen_tol
 */
void Grid::SetTol(const real_t refine_tol, const real_t coarsen_tol) {
    m_begin;
    //-------------------------------------------------------------------------
    rtol_ = refine_tol;
    ctol_ = coarsen_tol;
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief adapt = (refine or coarsen once) each block by one level, based on the given field
 * 
 * @param field the field used for the criterion
 */
void Grid::Adapt(Field* field) {
    m_begin;
    m_log("--> grid adaptation started... (interpolator: %s)", interp_->Identity().c_str());
    //-------------------------------------------------------------------------
    // store the criterion field
    tmp_ptr_ = reinterpret_cast<void*>(field);
    // compute the ghost needed by the interpolation of everyblock
    for (auto fid = fields_.begin(); fid != fields_.end(); fid++) {
        GhostPull(fid->second);
    }
    // delete the soon-to be outdated ghost and mesh
    DestroyGhost();
    // coarsen the needed block
    if (prof_ != nullptr) {
        prof_->Start("p4est_refcoarse");
    }
    // set the grid in the forest for the callback
    forest_->user_pointer = reinterpret_cast<void*>(this);
    p8est_coarsen_ext(forest_, 0, 0, cback_Interpolator, nullptr, cback_Interpolate);
    // refine the needed blocks
    p8est_refine_ext(forest_, 0, P8EST_MAXLEVEL, cback_Interpolator, nullptr, cback_Interpolate);
    // balance the partition
    if (prof_ != nullptr) {
        prof_->Stop("p4est_refcoarse");
        prof_->Start("p4est_balance");
    }
    p8est_balance_ext(forest_, P8EST_CONNECT_FULL, nullptr, cback_Interpolate);
    if (prof_ != nullptr) {
        prof_->Stop("p4est_balance");
        prof_->Start("p4est_partition_init");
    }
    // partition the grid
    Partitioner partition(&fields_, this, true);
    if (prof_ != nullptr) {
        prof_->Stop("p4est_partition_init");
        prof_->Start("p4est_partition_comm");
    }
    partition.Start(&fields_,M_FORWARD);
    partition.End(&fields_,M_FORWARD);
    if (prof_ != nullptr) {
        prof_->Stop("p4est_partition_comm");
    }
    // create a new ghost and mesh
    SetupGhost();
    // set the ghosting fields as non-valid
    for (auto fid = fields_.begin(); fid != fields_.end(); fid++) {
        fid->second->ghost_status(false);
    }
    // reset the forest pointer
    forest_->user_pointer = nullptr;
    //-------------------------------------------------------------------------
    m_log("--> grid adaptation done: now %ld blocks on %ld trees using %d ranks and %d threads", forest_->global_num_quadrants, forest_->trees->elem_count, forest_->mpisize, omp_get_max_threads());
    m_end;
}

/**
 * @brief adapt = (refine or coarsen once) each block recursivelly to reach the criterion imposed in the patch list.
 * 
 * @warning every existing field will be reset to 0 during its execution
 * 
 * A block is refined/coarsened if one part of it belongs to the patch. Simply touching the patch does not count as belonging (i.e.
 * the comparison is done using strict operators `<` and `>` ), see @ref cback_Patch().
 * 
 * @param patches the list of patch to match
 */
void Grid::Adapt(list<Patch>* patches) {
    m_begin;
    m_log("--> grid adaptation started... (using patches)");
    //-------------------------------------------------------------------------
    if (patches->size() == 0) {
        return;
    }
    // store the criterion patch list
    tmp_ptr_ = reinterpret_cast<void*>(patches);
    // no ghost is computed as no interpolation will be done
    // delete the soon-to be outdated ghost and mesh
    DestroyGhost();
    // set the grid in the forest for the callback
    forest_->user_pointer = reinterpret_cast<void*>(this);
    // we coarsen recursivelly
    m_verb("starting coarsening");
    if (prof_ != nullptr) {
        prof_->Start("p4est_refcoarse");
    }
    p8est_coarsen_ext(forest_, 1, 0, cback_Patch, nullptr, cback_AllocateOnly);
    // refine the needed blocks recursivelly
    p8est_refine_ext(forest_, 1, P8EST_MAXLEVEL, cback_Patch, nullptr, cback_AllocateOnly);
    // balance the partition
    if (prof_ != nullptr) {
        prof_->Stop("p4est_refcoarse");
        prof_->Start("p4est_balance");
    }
    p8est_balance_ext(forest_, P8EST_CONNECT_FULL, nullptr, cback_AllocateOnly);
    if (prof_ != nullptr) {
        prof_->Stop("p4est_balance");
        prof_->Start("p4est_partition_init");
    }
    // partition the grid
    Partitioner partition(&fields_, this, true);
    if (prof_ != nullptr) {
        prof_->Stop("p4est_partition_init");
        prof_->Start("p4est_partition_comm");
    }
    partition.Start(&fields_,M_FORWARD);
    partition.End(&fields_,M_FORWARD);
    if (prof_ != nullptr) {
        prof_->Stop("p4est_partition_comm");
    }
    // create a new ghost and mesh
    SetupGhost();
    // set the ghosting fields as non-valid
    for (auto fid = fields_.begin(); fid != fields_.end(); fid++) {
        fid->second->ghost_status(false);
    }
    // reset the forest pointer
    forest_->user_pointer = nullptr;
    //-------------------------------------------------------------------------
    m_log("--> grid adaptation done: now %ld blocks on %ld trees using %d ranks and %d threads", forest_->global_num_quadrants, forest_->trees->elem_count, forest_->mpisize, omp_get_max_threads());
    m_end;
}

/**
 * @brief iterates on the blocks and performs a simple @ref gbop_t operation using the forest structure and not the mesh
 * 
 * @warning for allocation and block management only. Use Operators (see operator.hpp) for computations
 * 
 * @param op 
 * @param field 
 */
void Grid::LoopOnGridBlock_(const gbop_t op, Field* field) const {
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
