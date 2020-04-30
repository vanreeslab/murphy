#include "grid.hpp"

#include <omp.h>
#include <p8est_extended.h>

#include "gridcallback.hpp"
#include "operator.hpp"
#include "partitioner.hpp"
#include "wavelet.hpp"

using std::string;
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
    // create a default interpolator
    interp_ = new Wavelet<5>();
    detail_ = new Wavelet<3>();

    // profiler
    prof_ = prof;

    // create the associated blocks
    p8est_iterate(forest_, NULL, NULL, cback_CreateBlock, NULL, NULL, NULL);

    // create the ghosts structure
    ghost_ = new Ghost(this,interp_);

    // init the profiler tracking
    if(prof_ != nullptr){
        prof->Create("ghost");
        prof->Create("mirror","ghost");
        prof->Create("send","ghost");
        prof->Create("recv","ghost");
        prof->Create("fill","ghost");
        prof->Create("adapt");
        prof->Create("get_ghost","adapt");
        prof->Create("coarsen","adapt");
        prof->Create("refine","adapt");
        prof->Create("balance","adapt");
        prof->Create("partition","adapt");
        prof->Create("re-setup","adapt");
    }
    //-------------------------------------------------------------------------
    m_log("uniform grid created with %ld blocks on %ld trees using %d ranks and %d threads", forest_->global_num_quadrants, forest_->trees->elem_count, forest_->mpisize, omp_get_max_threads());
    m_end;
}

/**
 * @brief Destroy the Grid, frees all the blocks and the fields contained (if not done)
 * 
 */
Grid::~Grid() {
    m_begin;
    //-------------------------------------------------------------------------
    // destroy the interpolator
    delete (interp_);
    delete (detail_);
    delete (ghost_);
    // destroy the remaining blocks
    p8est_iterate(forest_, NULL, NULL, cback_DestroyBlock, NULL, NULL, NULL);
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
 * @brief Pull ghost points (take the values from the neighbors): fill the mirror buffers and start the send MPI call
 * 
 * @warning this function is part of the advanced control feature
 * 
 * @param field the considered field
 * @param ida the dimension of the field which has to be send
 */
void Grid::GhostPullSend(Field* field, const sid_t ida) {
    m_begin;
    m_assert(0 <= ida && ida < field->lda(), "the ida is not within the field's limit");
    m_assert(field != nullptr, "the source field cannot be null");
    m_assert(IsAField(field), "the field does not belong to this grid");
    //-------------------------------------------------------------------------
    if (!field->ghost_status()) {
        if(prof_!=nullptr){
            prof_->Start("mirror");
        }
        ghost_->PushToMirror(field, ida);
        if(prof_!=nullptr){
            prof_->Stop("mirror");
            prof_->Start("send");
        }
        ghost_->MirrorToGhostSend();
        if(prof_!=nullptr){
            prof_->Stop("send");
        }
    }
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief Pull ghost points (take the values from the neighbors): receive the buffers.
 * After this function, the receive buffers are full and the send buffers are available for another send
 * 
 * @warning this function is part of the advanced control feature
 * 
 * @param field the considered field
 * @param ida the received dimension
 */
void Grid::GhostPullRecv(Field* field, const sid_t ida) {
    m_begin;
    m_assert(0 <= ida && ida < field->lda(), "the ida is not within the field's limit");
    m_assert(field != nullptr, "the source field cannot be null");
    m_assert(IsAField(field), "the field does not belong to this grid");
    //-------------------------------------------------------------------------
    if (!field->ghost_status()) {
        // receive the current communication, the mirrors are now free
        if(prof_!=nullptr){
            prof_->Start("recv");
        }
        ghost_->MirrorToGhostRecv();
        if(prof_!=nullptr){
            prof_->Stop("recv");
        }
    }
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief Pull ghost points (take the values from the neighbors): get the actual ghost points values from the received buffers
 * After this function, the reception buffers are available for another reception
 * 
 * @warning this function is part of the advanced control feature
 * 
 * @param field the considered field
 * @param ida the filled dimension
 */
void Grid::GhostPullFill(Field* field, const sid_t ida) {
    m_begin;
    m_assert(0 <= ida && ida < field->lda(), "the ida is not within the field's limit");
    m_assert(field != nullptr, "the source field cannot be null");
    m_assert(IsAField(field), "the field does not belong to this grid");
    m_assert(interp_ != nullptr, "the inteprolator cannot be null");
    //-------------------------------------------------------------------------
    if (!field->ghost_status()) {
        // receive the current communication, the mirrors are now free
        if(prof_!=nullptr){
            prof_->Start("fill");
        }
        ghost_->PullFromGhost(field, ida);
        if(prof_!=nullptr){
            prof_->Stop("fill");
        }
    }
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief Pull the ghost points (take the values from the neighbors)
 * 
 * It implements the overlapping between the send and the reception of every dimension in order
 * reduce the total computational time.
 * 
 * After this function, the ghost status of the field is set as up-to-date
 * 
 * @param field the field which requires the ghost
 */
void Grid::GhostPull(Field* field) {
    m_begin;
    m_assert(field != nullptr, "the source field cannot be null");
    //-------------------------------------------------------------------------
    // start the send in the first dimension
    GhostPullSend(field, 0);
    for (int ida = 1; ida < field->lda(); ida++) {
        // receive the previous dimension
        GhostPullRecv(field, ida - 1);
        // start the send for the next dimension
        GhostPullSend(field, ida);
        // fill the ghost values of the just-received information
        GhostPullFill(field, ida - 1);
    }
    GhostPullRecv(field, field->lda() - 1);
    // fill the ghost values of the just-received information
    GhostPullFill(field, field->lda() - 1);
    // set that everything is ready for the field
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
        delete (ghost_);
        ResetP4estGhostMesh();
        // set the grid in the forest for the callback
        forest_->user_pointer = reinterpret_cast<void*>(this);
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
        ghost_ = new Ghost(this,interp_);
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
        delete (ghost_);
        ResetP4estGhostMesh();
        // set the grid in the forest for the callback
        forest_->user_pointer = reinterpret_cast<void*>(this);
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
        ghost_ = new Ghost(this,interp_);
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
    m_log("--> grid adaptation started... (interpolator: %s)",detail_->Identity().c_str());
    //-------------------------------------------------------------------------
    // store the criterion field
    tmp_field_ = field;
    if(prof_ != nullptr){
        prof_->Start("adapt");
        prof_->Start("get_ghost");
    }
    // compute the ghost needed by the interpolation of everyblock
    for (auto fid = fields_.begin(); fid != fields_.end(); fid++) {
        GhostPull(fid->second);
    }
    if(prof_ != nullptr){
        prof_->Stop("get_ghost");
    }
    // delete the soon-to be outdated ghost and mesh
    delete (ghost_);
    ResetP4estGhostMesh();
    // set the grid in the forest for the callback
    forest_->user_pointer = reinterpret_cast<void*>(this);
    // coarsen the needed block
    if(prof_ != nullptr){
        prof_->Start("coarsen");
    }
    p8est_coarsen_ext(forest_, 0, 0, cback_Interpolator, nullptr, cback_Interpolate);
    if(prof_ != nullptr){
        prof_->Stop("coarsen");
        prof_->Start("refine");
    }
    // refine the needed blocks
    p8est_refine_ext(forest_, 0, P8EST_MAXLEVEL, cback_Interpolator, nullptr, cback_Interpolate);
    // balance the partition
    if(prof_ != nullptr){
        prof_->Stop("refine");
        prof_->Start("balance");
    }
    p8est_balance_ext(forest_, P8EST_CONNECT_FULL, NULL, cback_Interpolate);
    if(prof_ != nullptr){
        prof_->Stop("balance");
        prof_->Start("partition");
    }
    // partition the grid
    Partitioner partition(&fields_, this);
    partition.Start(&fields_);
    partition.End(&fields_);
    if(prof_ != nullptr){
        prof_->Stop("partition");
        prof_->Start("re-setup");
    }
    // create a new ghost and mesh
    SetupP4estGhostMesh();
    ghost_ = new Ghost(this,interp_);
    if(prof_ != nullptr){
        prof_->Stop("re-setup");
        prof_->Stop("adapt");
    }
    // set the ghosting fields as non-valid
    for (auto fid = fields_.begin(); fid != fields_.end(); fid++) {
        fid->second->ghost_status(false);
    }
    //-------------------------------------------------------------------------
    m_log("--> grid adaptation done: now %ld blocks on %ld trees using %d ranks and %d threads", forest_->global_num_quadrants, forest_->trees->elem_count, forest_->mpisize, omp_get_max_threads());
    m_end;
}

/**
 * @brief iterates on the blocks and performs a simple @ref bop_t operation using the forest structure and not the mesh
 * 
 * @warning for allocation and block management only. Use Operators (see operator.hpp) for computations
 * 
 * @param op 
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
