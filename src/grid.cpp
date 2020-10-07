#include "grid.hpp"

#include <omp.h>
#include <p8est_extended.h>

#include "gridcallback.hpp"
#include "partitioner.hpp"
#include "wavelet.hpp"

using std::list;
using std::string;

/**
 * @brief Construct a new Grid which is empty, only the interpolators have been associated
 * 
 */
Grid::Grid() : ForestGrid() {
    //-------------------------------------------------------------------------
    prof_  = nullptr;
    ghost_ = nullptr;
    // create a default interpolator
    interp_ = new Wavelet();
    //-------------------------------------------------------------------------
};

/**
 * @brief Construct a new Grid as a uniform grid, distributed among the cpus
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
    // m_profCreate(prof_, "p4est_refcoarse");
    // m_profCreate(prof_, "p4est balance");
    // m_profCreate(prof_, "partition init");
    // m_profCreate(prof_, "partition comm");
    // m_profCreateParent(prof_, "p4est_refcoarse", "cback_interpolate");

    // create a default interpolator
    interp_ = new Wavelet();
    // create the associated blocks
    p8est_iterate(p4est_forest_, NULL, NULL, cback_CreateBlock, NULL, NULL, NULL);
    // partition the grid to have compatible grid
    Partitioner part = Partitioner(&fields_, this, true);
    part.Start(&fields_, M_FORWARD);
    part.End(&fields_, M_FORWARD);
    // setup the ghost stuctures as the mesh will not change anymore
    m_verb("The grid is partitioned -> let's build the ghost now");
    SetupGhost();
    //-------------------------------------------------------------------------
    m_log("uniform grid created with %ld blocks on %ld trees using %d ranks and %d threads", p4est_forest_->global_num_quadrants, p4est_forest_->trees->elem_count, p4est_forest_->mpisize, omp_get_max_threads());
    m_end;
}

/**
 * @brief Copy the ForestGrid part from a grid
 * 
 * @param grid the source grid
 */
void Grid::CopyFrom(Grid* grid) {
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
        m_verb("dealloc the interp");
        delete (interp_);
    }
    // destroy the ghosts, they are mine as well
    DestroyGhost();
    // destroy the remaining blocks
    if (is_connect_owned_) {
        p8est_iterate(p4est_forest_, NULL, NULL, cback_DestroyBlock, NULL, NULL, NULL);
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
    m_verb("starting the Ghost construction");
    ghost_ = new Ghost(this, interp_, prof_);
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
        m_verb("dealloc the ghost");
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
        memsize += p4est_forest_->local_num_quadrants * (M_N * M_N * M_N) * fid->second->lda() * sizeof(real_t);
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
    return p4est_forest_->local_num_quadrants * (M_N * M_N * M_N);
}

/**
 * @brief returns the number of global number of degree of freedom for one field containing one dimension
 * 
 * @return size_t 
 */
size_t Grid::GlobalNumDof() const {
    return p4est_forest_->global_num_quadrants * (M_N * M_N * M_N);
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
        // LoopOnGridBlock_(&GridBlock::AddField, field);
        // DoOpTree(&CallGridBlockMemFuncField, this, field, &GridBlock::AddField);
        DoOpTree(nullptr, &GridBlock::AddField, this, field);
        m_log("field %s has been added to the grid", field->name().c_str());
    } else {
        m_log("field %s is already in the grid", field->name().c_str());
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
        // remove the field to everyblock
        // DoOpTree(&CallGridBlockMemFuncField, this, field, &GridBlock::DeleteField);
        DoOpTree(nullptr, &GridBlock::DeleteField, this, field);
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
        m_assert(iter->first == iter->second->name(), "the key of the map must be the name of the field");
    }
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief Pull the ghost points (take values from neighbors to fill my ghost points) - start the comm
 * 
 * @param field the considered field
 * @param ida the dimension of the field which has to be send
 */
void Grid::GhostPull_Post(Field* field, const sid_t ida) {
    m_begin;
    m_assert(0 <= ida && ida < field->lda(), "the ida is not within the field's limit");
    m_assert(field != nullptr, "the source field cannot be null");
    m_assert(IsAField(field), "the field does not belong to this grid");
    m_assert(ghost_ != nullptr, "The ghost structure is not valid, unable to use it");
    //-------------------------------------------------------------------------
    if (!field->ghost_status()) {
        ghost_->PullGhost_Post(field, ida);
    }
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief Pull the ghost points (take values from neighbors to fill my ghost points) - end the comm
 * 
 * @param field the considered field
 * @param ida the received dimension
 */
void Grid::GhostPull_Wait(Field* field, const sid_t ida) {
    m_begin;
    m_assert(0 <= ida && ida < field->lda(), "the ida is not within the field's limit");
    m_assert(field != nullptr, "the source field cannot be null");
    m_assert(IsAField(field), "the field does not belong to this grid");
    m_assert(ghost_ != nullptr, "The ghost structure is not valid, unable to use it");
    //-------------------------------------------------------------------------
    if (!field->ghost_status()) {
        ghost_->PullGhost_Wait(field, ida);
    }
    // put the field as correct
    field->ghost_status(true);
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
    m_profStart(prof_, "pull ghost");
    for (int ida = 0; ida < field->lda(); ida++) {
        GhostPull_Post(field, ida);
        GhostPull_Wait(field, ida);
    }
    m_profStop(prof_, "pull ghost");
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
    m_profStart(prof_,"refine");
    // init the profiler
    m_profInit(prof_,"p4est refine");
    m_profInitLeave(prof_,"wavelet interpolation");
    m_profLeave(prof_,"p4est refine");
    m_profInit(prof_,"p4est balance");
    m_profInitLeave(prof_,"wavelet interpolation");
    m_profLeave(prof_, "p4est balance");
    //................................................

    // set the grid in the forest for the callback
    p4est_forest_->user_pointer = reinterpret_cast<void*>(this);

    // we create the new blocks
    for (int id = 0; id < delta_level; id++) {
        // unlock all the blocks
        // LoopOnGridBlock_(&GridBlock::unlock, nullptr);
        // DoOpTree(&CallGridBlockMemFuncEmpty, this, &GridBlock::unlock);
        DoOpTree(nullptr, &GridBlock::unlock, this);

        // compute the ghost needed by the interpolation
        for (auto fid = fields_.begin(); fid != fields_.end(); fid++) {
            GhostPull(fid->second);
        }
        // delete the soon-to be outdated ghost and mesh
        DestroyGhost();

        // do the p4est interpolation by callback
        m_profStart(prof_, "p4est refine");
        p8est_refine_ext(p4est_forest_, 0, P8EST_MAXLEVEL, cback_Yes, nullptr, cback_Interpolate);
        m_profStop(prof_, "p4est refine");

        // balance the partition
        m_profStart(prof_, "p4est balance");
        p8est_balance_ext(p4est_forest_, P8EST_CONNECT_FULL, NULL, cback_Interpolate);
        m_profStop(prof_, "p4est balance");

        // partition the grid
        m_profStart(prof_, "partition init");
        Partitioner partition(&fields_, this, true);
        m_profStop(prof_, "partition init");
        m_profStart(prof_, "partition comm");
        partition.Start(&fields_, M_FORWARD);
        partition.End(&fields_, M_FORWARD);
        m_profStop(prof_, "partition comm");

        // create a new ghost and mesh
        SetupGhost();
        // set the ghosting fields as non-valid
        for (auto fid = fields_.begin(); fid != fields_.end(); fid++) {
            fid->second->ghost_status(false);
        }
    }

    p4est_forest_->user_pointer = nullptr;
    m_profStop(prof_,"refine");
    //-------------------------------------------------------------------------
    m_assert(p4est_forest_->user_pointer == nullptr, "we must reset the user_pointer to null");
    m_log("refined grid created with %ld blocks on %ld trees using %d ranks and %d threads", p4est_forest_->global_num_quadrants, p4est_forest_->trees->elem_count, p4est_forest_->mpisize, omp_get_max_threads());
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
    m_profStart(prof_,"coarsen");
    //................................................
    m_profInit(prof_,"p4est coarsen");
    m_profInitLeave(prof_,"wavelet interpolation");
    m_profLeave(prof_,"p4est coarsen");
    m_profInit(prof_,"p4est balance");
    m_profInitLeave(prof_,"wavelet interpolation");
    m_profLeave(prof_, "p4est balance");
    //................................................

    // set the grid in the forest for the callback
    p4est_forest_->user_pointer = reinterpret_cast<void*>(this);

    // we create the new blocks
    for (int id = 0; id < delta_level; id++) {
        // unlock all the blocks
        // LoopOnGridBlock_(&GridBlock::unlock, nullptr);
        // DoOpTree(&CallGridBlockMemFuncEmpty, this, &GridBlock::unlock);
        DoOpTree(nullptr, &GridBlock::unlock, this);

        // compute the ghost needed by the interpolation
        for (auto fid = fields_.begin(); fid != fields_.end(); fid++) {
            // set the working field before entering the callback
            // working_callback_field_ = fid->second;
            // get the ghosts needed by the interpolation
            GhostPull(fid->second);
        }
        // delete the soon-to be outdated ghost and mesh
        DestroyGhost();

        // do the p4est interpolation by callback
        m_profStart(prof_, "p4est coarsen");
        p8est_coarsen_ext(p4est_forest_, 0, 0, cback_Yes, nullptr, cback_Interpolate);
        m_profStop(prof_, "p4est coarsen");

        // balance the partition
        m_profStart(prof_, "p4est balance");
        p8est_balance_ext(p4est_forest_, P8EST_CONNECT_FULL, NULL, cback_Interpolate);
        m_profStop(prof_, "p4est balance");

        // partition the grid
        m_profStart(prof_, "partition init");
        Partitioner partition(&fields_, this, true);
        m_profStop(prof_, "partition init");
        m_profStart(prof_, "partition comm");
        partition.Start(&fields_, M_FORWARD);
        partition.End(&fields_, M_FORWARD);
        m_profStop(prof_, "partition comm");

        // create a new ghost and mesh
        SetupGhost();
        // set the ghosting fields as non-valid
        for (auto fid = fields_.begin(); fid != fields_.end(); fid++) {
            fid->second->ghost_status(false);
        }
    }

    p4est_forest_->user_pointer = nullptr;
    m_profStop(prof_,"coarsen");
    //-------------------------------------------------------------------------
    m_assert(p4est_forest_->user_pointer == nullptr, "we must reset the user_pointer to null");
    m_log("coarsened grid created with %ld blocks on %ld trees using %d ranks and %d threads", p4est_forest_->global_num_quadrants, p4est_forest_->trees->elem_count, p4est_forest_->mpisize, omp_get_max_threads());
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
    m_assert(refine_tol > coarsen_tol, "the refinement tolerance must be > the coarsening tolerance");
    //-------------------------------------------------------------------------
    rtol_ = refine_tol;
    ctol_ = coarsen_tol;
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief adapt = (refine or coarsen !once!) each block by one level, based on the given field
 * 
 * @param field the field used for the criterion
 */
void Grid::Adapt(Field* field) {
    m_begin;
    m_assert(IsAField(field), "the field must already exist on the grid!");
    m_assert(cback_criterion_field_ == nullptr, "the pointer `cback_criterion_field` must be  null");
    m_assert(cback_interpolate_ptr_ == nullptr, "the pointer `cback_interpolate_ptr` must be  null");
    m_assert(p4est_forest_->user_pointer == nullptr, "we must start with the user_pointer to null");
    m_assert(!recursive_adapt(), "we cannot refine recursivelly here");
    m_log("--> grid adaptation started... (interpolator: %s)", interp_->Identity().c_str());
    //---------------------------------------------s----------------------------
    m_profStart(prof_, "adapt field");
    //................................................
    m_profInit(prof_, "p4est coarsen");
    m_profInitLeave(prof_, "wavelet interpolation");
    m_profInitLeave(prof_, "wavelet criterion");
    m_profLeave(prof_, "p4est coarsen");
    m_profInit(prof_, "p4est refine");
    m_profInitLeave(prof_, "wavelet interpolation");
    m_profInitLeave(prof_, "wavelet criterion");
    m_profLeave(prof_, "p4est refine");
    m_profInit(prof_, "p4est balance");
    m_profInitLeave(prof_, "wavelet interpolation");
    m_profLeave(prof_, "p4est balance");
    //................................................

    // unlock all the blocks
    DoOpTree(nullptr, &GridBlock::unlock, this);

    // compute the ghost needed by the interpolation of every other field in the grid
    for (auto fid = fields_.begin(); fid != fields_.end(); fid++) {
        GhostPull(fid->second);
    }
    // delete the soon-to be outdated ghost and mesh
    DestroyGhost();

    // store the criterion field
    cback_criterion_field_      = reinterpret_cast<void*>(field);
    p4est_forest_->user_pointer = reinterpret_cast<void*>(this);

    // coarsen the needed block
    m_profStart(prof_, "p4est coarsen");
    p8est_coarsen_ext(p4est_forest_, recursive_adapt(), 0, cback_Interpolator, nullptr, cback_Interpolate);
    m_profStop(prof_, "p4est coarsen");
    m_profStart(prof_, "p4est refine");
    p8est_refine_ext(p4est_forest_, recursive_adapt(), P8EST_MAXLEVEL, cback_Interpolator, nullptr, cback_Interpolate);
    m_profStop(prof_, "p4est refine");
    // balance the partition
    m_profStart(prof_, "p4est balance");
    p8est_balance_ext(p4est_forest_, P8EST_CONNECT_FULL, nullptr, cback_Interpolate);
    m_profStop(prof_, "p4est balance");

    // partition the grid
    m_profStart(prof_, "partition init");
    Partitioner partition(&fields_, this, true);
    m_profStop(prof_, "partition init");
    m_profStart(prof_, "partition comm");
    partition.Start(&fields_, M_FORWARD);
    partition.End(&fields_, M_FORWARD);
    m_profStop(prof_, "partition comm");

    // create a new ghost and mesh
    SetupGhost();

    // set the ghosting fields as non-valid
    for (auto fid = fields_.begin(); fid != fields_.end(); fid++) {
        fid->second->ghost_status(false);
    }
    // reset the forest pointer
    cback_criterion_field_      = nullptr;
    p4est_forest_->user_pointer = nullptr;
    m_profStop(prof_, "adapt field");
    //-------------------------------------------------------------------------
    m_assert(cback_criterion_field_ == nullptr, "the pointer `cback_criterion_field_` must be  null");
    m_assert(cback_interpolate_ptr_ == nullptr, "the pointer `cback_interpolate_ptr` must be  null");
    m_assert(p4est_forest_->user_pointer == nullptr, "we must reset the user_pointer to null");
    m_log("--> grid adaptation done: now %ld blocks on %ld trees using %d ranks and %d threads", p4est_forest_->global_num_quadrants, p4est_forest_->trees->elem_count, p4est_forest_->mpisize, omp_get_max_threads());
    m_end;
}
/**
 * @brief Adapt the grid given an analytical boundary condition for the only field present in the grid
 * 
 * @warning only field will be set to the value, the rest will be set to 0.0 or kept as such
 * 
 * @param expression 
 */
void Grid::Adapt(Field* field, SetValue* expression) {
    m_begin;
    m_assert(IsAField(field), "the field must already exist on the grid!");
    m_assert(cback_criterion_field_ == nullptr, "the pointer `cback_criterion_field_` must be  null");
    m_assert(cback_interpolate_ptr_ == nullptr, "the pointer `cback_interpolate_ptr` must be  null");
    m_assert(p4est_forest_->user_pointer == nullptr, "we must start with the user_pointer to null");
    m_log("--> grid adaptation started... (using analytical expression and recursive = %d)", recursive_adapt());
    //-------------------------------------------------------------------------
    m_profStart(prof_, "adapt expression");
    //................................................
    m_profInit(prof_, "p4est coarsen");
    m_profInitLeave(prof_, "wavelet criterion");
    m_profInitLeave(prof_, "fill values");
    m_profLeave(prof_, "p4est coarsen");
    m_profInit(prof_, "p4est refine");
    m_profInitLeave(prof_, "wavelet criterion");
    m_profInitLeave(prof_, "fill values");
    m_profLeave(prof_, "p4est refine");
    m_profInit(prof_, "p4est balance");
    m_profInitLeave(prof_, "fill values");
    m_profLeave(prof_, "p4est balance");
    //................................................
    // unlock all the blocks
    DoOpTree(nullptr, &GridBlock::unlock, this);

    // delete the soon-to be outdated ghost and mesh
    DestroyGhost();

    // store the operator and the grid for access in the callback function
    cback_criterion_field_      = reinterpret_cast<void*>(field);
    cback_interpolate_ptr_      = reinterpret_cast<void*>(expression);
    p4est_forest_->user_pointer = reinterpret_cast<void*>(this);
    m_verb("Storing the field %s for refinement", field->name().c_str());

    // apply the operator to get the starting value
    (*expression)(this, field);
    m_assert(field->ghost_status(), "the ghost status should be valid here...");

    // coarsen + refine the needed blocks recursivelly using the wavelet criterion and the operator fill
    m_profStart(prof_, "p4est coarsen");
    p8est_coarsen_ext(p4est_forest_, recursive_adapt(), 0, cback_Interpolator, nullptr, cback_ValueFill);
    m_profStop(prof_, "p4est coarsen");
    m_profStart(prof_, "p4est refine");
    p8est_refine_ext(p4est_forest_, recursive_adapt(), P8EST_MAXLEVEL, cback_Interpolator, nullptr, cback_ValueFill);
    m_profStop(prof_, "p4est refine");

    // balance the partition
    m_profStart(prof_, "p4est balance");
    p8est_balance_ext(p4est_forest_, P8EST_CONNECT_FULL, nullptr, cback_ValueFill);
    m_profStop(prof_, "p4est balance");

    // partition the grid
    m_profStart(prof_, "partition init");
    Partitioner partition(&fields_, this, true);
    m_profStop(prof_, "partition init");
    m_profStart(prof_, "partition comm");
    partition.Start(&fields_, M_FORWARD);
    partition.End(&fields_, M_FORWARD);
    m_profStop(prof_, "partition comm");

    // create a new ghost and mesh
    SetupGhost();
    // set the ghosting fields as non-valid
    for (auto fid = fields_.begin(); fid != fields_.end(); fid++) {
        fid->second->ghost_status(false);
    }
    // we modified one block after another, so we set the ghost value from the SetValue
    field->ghost_status(expression->do_ghost());

    // reset the forest pointer
    cback_criterion_field_      = nullptr;
    cback_interpolate_ptr_      = nullptr;
    p4est_forest_->user_pointer = nullptr;

    // reset the recursive to false
    SetRecursiveAdapt(false);
    m_profStop(prof_, "adapt expression");
    //-------------------------------------------------------------------------
    m_assert(cback_criterion_field_ == nullptr, "the pointer `cback_criterion_field_` must be  null");
    m_assert(cback_interpolate_ptr_ == nullptr, "the pointer `cback_interpolate_ptr` must be  null");
    m_assert(p4est_forest_->user_pointer == nullptr, "we must reset the user_pointer to null");
    m_log("--> grid adaptation done: now %ld blocks on %ld trees using %d ranks and %d threads", p4est_forest_->global_num_quadrants, p4est_forest_->trees->elem_count, p4est_forest_->mpisize, omp_get_max_threads());
    m_end;
}

/**
 * @brief adapt = (refine or coarsen) each block recursivelly to reach the criterion imposed in the patch list.
 * 
 * @warning no interpolation will be performed, i.e. every existing field will be reset to 0 during its execution
 * 
 * A block is refined/coarsened if one part of it belongs to the patch. Simply touching the patch does not count as belonging (i.e.
 * the comparison is done using strict operators `<` and `>` ), see @ref cback_Patch().
 * 
 * @param patches the list of patch to match
 */
void Grid::Adapt(list<Patch>* patches) {
    m_begin;
    m_assert(cback_criterion_field_ == nullptr, "the pointer `cback_criterion_field_` must be  null");
    m_assert(cback_interpolate_ptr_ == nullptr, "the pointer `cback_interpolate_ptr` must be  null");
    m_assert(p4est_forest_->user_pointer == nullptr, "we must start with the user_pointer to null");
    m_log("--> grid adaptation started... (using patches and recursive = %d)", recursive_adapt());
    //-------------------------------------------------------------------------
    m_profStart(prof_, "adapt patch");
    //................................................
    m_profInit(prof_, "p4est coarsen");
    m_profInitLeave(prof_, "allocate only");
    m_profLeave(prof_, "p4est coarsen");
    m_profInit(prof_, "p4est refine");
    m_profInitLeave(prof_, "allocate only");
    m_profLeave(prof_, "p4est refine");
    m_profInit(prof_, "p4est balance");
    m_profInitLeave(prof_, "allocate only");
    m_profLeave(prof_, "p4est balance");
    //................................................

    if (patches->size() == 0) {
        return;
    }
    // unlock all the blocks
    // LoopOnGridBlock_(&GridBlock::unlock, nullptr);
    // DoOpTree(&CallGridBlockMemFuncEmpty, this, &GridBlock::unlock);
    DoOpTree(nullptr, &GridBlock::unlock, this);

    // delete the soon-to be outdated ghost and mesh
    DestroyGhost();

    // store the criterion patch list
    cback_criterion_field_      = reinterpret_cast<void*>(patches);
    p4est_forest_->user_pointer = reinterpret_cast<void*>(this);

    // coarsen + refine the needed blocks recursivelly
    m_profStart(prof_, "p4est coarsen");
    p8est_coarsen_ext(p4est_forest_, recursive_adapt(), 0, cback_Patch, nullptr, cback_AllocateOnly);
    m_profStop(prof_, "p4est coarsen");
    m_profStart(prof_, "p4est refine");
    p8est_refine_ext(p4est_forest_, recursive_adapt(), P8EST_MAXLEVEL, cback_Patch, nullptr, cback_AllocateOnly);
    m_profStop(prof_, "p4est refine");

    // balance the partition
    m_profStart(prof_, "p4est balance");
    p8est_balance_ext(p4est_forest_, P8EST_CONNECT_FULL, nullptr, cback_AllocateOnly);
    m_profStop(prof_, "p4est balance");

    // partition the grid
    m_profStart(prof_, "partition init");
    Partitioner partition(&fields_, this, true);
    m_profStop(prof_, "partition init");
    m_profStart(prof_, "partition comm");
    partition.Start(&fields_, M_FORWARD);
    partition.End(&fields_, M_FORWARD);
    m_profStop(prof_, "partition comm");

    // create a new ghost and mesh
    SetupGhost();
    // set the ghosting fields as non-valid
    for (auto fid = fields_.begin(); fid != fields_.end(); fid++) {
        fid->second->ghost_status(false);
    }
    // reset the forest pointer
    cback_criterion_field_      = nullptr;
    p4est_forest_->user_pointer = nullptr;

    // reset the recursive to false
    SetRecursiveAdapt(false);

    m_profStop(prof_, "adapt patch");
    //-------------------------------------------------------------------------
    m_assert(cback_criterion_field_ == nullptr, "the pointer `cback_criterion_field` must be  null");
    m_assert(cback_interpolate_ptr_ == nullptr, "the pointer `cback_interpolate_ptr` must be  null");
    m_assert(p4est_forest_->user_pointer == nullptr, "we must reset the user_pointer to null");
    m_log("--> grid adaptation done: now %ld blocks on %ld trees using %d ranks and %d threads", p4est_forest_->global_num_quadrants, p4est_forest_->trees->elem_count, p4est_forest_->mpisize, omp_get_max_threads());
    m_end;
}
