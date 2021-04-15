#include "grid.hpp"

#include <omp.h>
#include <p8est_extended.h>

#include "partitioner.hpp"
#include "wavelet/interpolating_wavelet.hpp"

using std::list;
using std::string;

/**
 * @brief Construct a new Grid which is empty, only the Wavelets have been associated
 * 
 */
Grid::Grid() : ForestGrid(), prof_(nullptr), ghost_(nullptr), interp_(nullptr){};

/**
 * @brief Construct a new Grid as a uniform grid, distributed among the cpus
 * 
 * Initialize the ForestGrid, the Wavelet, the profiler, the ghosts and the grid blocks.
 * The grid is defined as a set of l[0]xl[1]xl[2] octrees, each of them refined up to level ilvl.
 * 
 * @param ilvl the starting initialization level
 * @param isper defines the periodicity of the dimension (cannot be changed afterwards)
 * @param l the aspect ratio of the domain (integers), defining how many trees are used in each direction
 * @param comm the MPI communicator used
 * @param prof the profiler pointer if any (can be nullptr)
 */
Grid::Grid(const level_t ilvl, const bool isper[3], const lid_t l[3], MPI_Comm comm, const m_ptr<Prof>& prof)
    : ForestGrid(ilvl, isper, l, sizeof(GridBlock*), comm) {
    m_begin;
    //-------------------------------------------------------------------------
    // profiler
    prof_ = prof();
    // create a default Wavelet -> default is M_WAVELET_N and M_WAVELET_NT
    interp_ = new InterpolatingWavelet();

    // create the associated blocks
    p8est_iterate(p4est_forest_, NULL, NULL, cback_CreateBlock, NULL, NULL, NULL);

    // partition the grid to have compatible grid
    Partitioner part = Partitioner(&fields_, this, true);
    part.Start(&fields_, M_FORWARD);
    part.End(&fields_, M_FORWARD);

    // setup the ghost stuctures as the mesh will not change anymore
    SetupMeshGhost();
    SetupAdapt();
    //-------------------------------------------------------------------------
    m_log("uniform grid created with %ld blocks on %ld trees using %d ranks and %d threads", p4est_forest_->global_num_quadrants, p4est_forest_->trees->elem_count, p4est_forest_->mpisize, omp_get_max_threads());
    m_end;
}

/**
 * @brief Copy the ForestGrid part from a grid
 * 
 * @param grid the source grid
 */
void Grid::CopyFrom(m_ptr<const Grid> grid) {
    m_begin;
    //-------------------------------------------------------------------------
    this->ForestGrid::CopyFrom(grid);
    // copy the field mapping
    for (auto iter = grid->FieldBegin(); iter != grid->FieldEnd(); iter++) {
        string name   = iter->first;
        Field* fid    = iter->second();
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
    // destroy the Wavelet and the details they are mine
    if (interp_ != nullptr) {
        m_verb("dealloc the interp");
        delete (interp_);
    }
    // destroy the ghosts, they are mine as well
    DestroyMeshGhost();
    DestroyAdapt();
    // destroy the remaining blocks
    if (is_connect_owned_) {
        p8est_iterate(p4est_forest_, NULL, NULL, cback_DestroyBlock, NULL, NULL, NULL);
    }
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief update the min and max level taken into account for the next refinement (does not change the current state of the grid to match the new requirements!)
 * 
 * the function fails if the existing levels on the grid do not fit into the new bounds (it avoid blocks blocked at a given level).
 * 
 */
void Grid::level_limit(const level_t min, const level_t max) {
    m_assert(min <= max, "the levels must be %d <= %d", min, max);
    m_assert(this->MinLevel() >= min, "trying to impose a min level = %d while blocks exist on a lower one = %d",min,this->MinLevel());
    m_assert(this->MaxLevel() <= max, "trying to impose a max level = %d while blocks exist on a higher one = %d",max,this->MaxLevel());
    //-------------------------------------------------------------------------
    level_limit_min_ = min;
    level_limit_max_ = max;
    m_verb("limit leves are now %d to %d ", min, max);
    //-------------------------------------------------------------------------
};

/**
 * @brief setup the Ghost structure when the mesh is not going to change anymore
 * 
 * @warning this function cannot be called on on existing structure
 * 
 */
void Grid::SetupMeshGhost() {
    m_begin;
    m_assert(ghost_ == nullptr, "cannot create something that already exists");
    //-------------------------------------------------------------------------
    // create the forestGrid part
    this->SetupP4estMeshAndGhost();
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
void Grid::DestroyMeshGhost() {
    m_begin;
    //-------------------------------------------------------------------------
    if (ghost_ != nullptr) {
        m_verb("dealloc the ghost");
        delete (ghost_);
        ghost_ = nullptr;
    }
    this->DestroyP4estMeshAndGhost();
    //-------------------------------------------------------------------------
    m_end;
}

void Grid::SetupAdapt() {
    m_begin;
    m_assert(is_mesh_valid(), "the mesh must be valid");
    //-------------------------------------------------------------------------
    // MPI_Info info;
    // MPI_Info_create(&info);
    // MPI_Info_set(info, "no_locks", "true");

    // MPI_Aint status_mem_size = local_num_quadrants() * sizeof(short);
    // MPI_Win_create(neighbor_status_window_, status_mem_size, sizeof(short), info, mpicomm(), &neighbor_status_window_);
    // neighbor_status_ = reinterpret_cast<short*>(m_calloc(status_mem_size));

    // MPI_Info_free(&info);
    //-------------------------------------------------------------------------
    m_end;
}

void Grid::DestroyAdapt() {
    m_begin;
    // //-------------------------------------------------------------------------
    // if (coarsen_status_ != nullptr) {
    //     m_free(coarsen_status_);
    //     coarsen_status_ = nullptr;
    // }
    // if (coarsen_status_window_ != MPI_WIN_NULL) {
    //     MPI_Win_free(&coarsen_status_window_);
    //     coarsen_status_window_ = MPI_WIN_NULL;
    // }
    //-------------------------------------------------------------------------
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
    for (auto fid = fields_.cbegin(); fid != fields_.cend(); fid++) {
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

bool Grid::IsAField(const m_ptr<const Field> field) const {
    std::string key = field->name();
    return (fields_.find(key) != fields_.end());
}

/**
 * @brief register a field in the current grid
 *
 * It loops on the blocks to allocate the needed memory
 * 
 * @param field 
 */
void Grid::AddField(m_ptr<Field> field) {
    m_begin;
    m_assert(!field.IsOwned(), "The field cannot be owned as it has not been created here");
    //-------------------------------------------------------------------------
    if (!IsAField(field)) {
        // add the field
        fields_[field->name()] = field;
        // allocate the field on every block
        DoOpTree(nullptr, &GridBlock::AddField, this, field);
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
void Grid::DeleteField(m_ptr<const Field> field) {
    m_begin;
    //-------------------------------------------------------------------------
    if (IsAField(field)) {
        // remove the field to everyblock
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
void Grid::ResetFields(m_ptr<const std::map<string, m_ptr<Field> > > fields){
    m_begin;
    //-------------------------------------------------------------------------
    // clear the current map
    fields_.clear();
    // copy the new one
    for (auto iter = fields->cbegin(); iter != fields->cend(); iter++) {
        fields_[iter->first] = iter->second();

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
void Grid::GhostPull_Post(m_ptr<const Field> field, const sid_t ida) const {
    m_begin;
    m_assert(0 <= ida && ida < field->lda(), "the ida is not within the field's limit");
    m_assert(!field.IsEmpty(), "the source field cannot be empty");
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
void Grid::GhostPull_Wait(m_ptr<const Field> field, const sid_t ida) const {
    m_begin;
    m_assert(0 <= ida && ida < field->lda(), "the ida is not within the field's limit");
    m_assert(!field.IsEmpty(), "the source field cannot be empty");
    m_assert(IsAField(field), "the field does not belong to this grid");
    m_assert(ghost_ != nullptr, "The ghost structure is not valid, unable to use it");
    //-------------------------------------------------------------------------
    if (!field->ghost_status()) {
        ghost_->PullGhost_Wait(field, ida);
    }
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief Pull the ghost points (take the values from the neighbors): after this function, the ghost status of the field is set as up-to-date
 * 
 * @param field the field which requires the ghost
 */
void Grid::GhostPull(m_ptr<Field> field) const {
    m_begin;
    m_assert(!field.IsEmpty(), "the source field cannot be empty");
    m_assert(ghost_ != nullptr, "The ghost structure is not valid, unable to use it");
    //-------------------------------------------------------------------------
    // start the send in the first dimension
    m_verb("ghost check: field <%s> is %s", field->name().c_str(), field->ghost_status() ? "OK" : "to be computed");
    m_profStart(prof_, "pull ghost");
    for (lda_t ida = 0; ida < field->lda(); ++ida) {
        m_verb("ghosting post field <%s> in dir %d", field->name().c_str(), ida);
        GhostPull_Post(field, ida);
        m_verb("ghosting wait field <%s> in dir %d", field->name().c_str(), ida);
        GhostPull_Wait(field, ida);
    }
    m_profStop(prof_, "pull ghost");
    // // set that everything is ready for the field
    field->ghost_status(true);
    //-------------------------------------------------------------------------
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
 * @brief Refine the grid if the criterion says so on the field
 * 
 * Interpolate the field (which are not temporary) to match the new mesh structure
 * 
 * @param field the field to use for the criterion computation
 */
void Grid::Refine(m_ptr<Field> field) {
    m_begin;
    m_assert(IsAField(field), "the field must already exist on the grid!");
    m_assert(!recursive_adapt(), "we cannot refine recursivelly here");
    //-------------------------------------------------------------------------
    // compute the ghost needed by the interpolation of every other field in the grid
    for (auto fid = fields_.begin(); fid != fields_.end(); fid++) {
        m_ptr<Field> cur_field = fid->second;
        if (!cur_field->is_temp()) {
            GhostPull(cur_field);
        }
    }

    // Adapt(reinterpret_cast<void*>(field), nullptr, nullptr, &cback_WaveDetail, &cback_Interpolate);
    // Adapt(reinterpret_cast<void*>(field), nullptr, nullptr, &cback_WaveDetail, &cback_UpdateDependency);
    // Adapt(field, nullptr, &cback_StatusCheck, reinterpret_cast<void*>(field()), cback_UpdateDependency, nullptr);
    AdaptMagic(field, nullptr, nullptr, &cback_StatusCheck, reinterpret_cast<void*>(field()), cback_UpdateDependency, nullptr);

    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief Coarsen the grid if the criterion says so on the field
 * 
 * Interpolate the field (which are not temporary) to match the new mesh structure
 * 
 * @param field the field to use for the criterion computation
 */
void Grid::Coarsen(m_ptr<Field> field) {
    m_begin;
    m_assert(IsAField(field), "the field must already exist on the grid!");
    m_assert(!recursive_adapt(), "we cannot refine recursivelly here");
    //-------------------------------------------------------------------------
    // compute the ghost needed by the interpolation of every other field in the grid
    for (auto fid = fields_.begin(); fid != fields_.end(); fid++) {
        m_ptr<Field> cur_field = fid->second;
        if (!cur_field->is_temp()) {
            GhostPull(cur_field);
        }
    }

    // Adapt(reinterpret_cast<void*>(field), nullptr, &cback_WaveDetail, nullptr, &cback_Interpolate);
    // Adapt(reinterpret_cast<void*>(field), nullptr, &cback_WaveDetail, nullptr, &cback_UpdateDependency);
    // Adapt(field, &cback_StatusCheck, nullptr, reinterpret_cast<void*>(field()), cback_UpdateDependency, nullptr);
    AdaptMagic(field, nullptr, &cback_StatusCheck, nullptr, reinterpret_cast<void*>(field()), cback_UpdateDependency, nullptr);

    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief adapt each block by one level, based on the given field and on the detail coefficient
 * 
 * Interpolate the field (which are not temporary) to match the new mesh structure
 * 
 * @param field the field used for the criterion
 */
void Grid::Adapt(m_ptr<Field> field) {
    m_begin;
    m_assert(IsAField(field), "the field must already exist on the grid!");
    m_assert(!recursive_adapt(), "we cannot refine recursivelly here");
    //-------------------------------------------------------------------------
    // compute the ghost needed by the interpolation of every other field in the grid
    for (auto fid = fields_.begin(); fid != fields_.end(); fid++) {
        m_ptr<Field> cur_field = fid->second;
        if (!cur_field->is_temp()) {
            GhostPull(cur_field);
        }
    }

    // Adapt(reinterpret_cast<void*>(field), nullptr, &cback_WaveDetail, &cback_WaveDetail, &cback_Interpolate);
    // Adapt(reinterpret_cast<void*>(field), nullptr, &cback_WaveDetail, &cback_WaveDetail, &cback_UpdateDependency);
    AdaptMagic(field, nullptr, &cback_StatusCheck, &cback_StatusCheck, reinterpret_cast<void*>(field()), &cback_UpdateDependency, nullptr);

    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief Adapt the grid given an analytical expression for the designated field
 * 
 * @warning For any other field than the one given in argument, no value is guaranteed
 * 
 * @param field 
 * @param expression 
 */
void Grid::Adapt(m_ptr<Field> field, m_ptr<SetValue> expression) {
    m_begin;
    m_assert(IsAField(field), "the field must already exist on the grid!");
    m_assert(!(recursive_adapt() && !expression->do_ghost()), "in recursive mode, I need to fill the ghost values while using the analytical expression");
    //-------------------------------------------------------------------------
    // apply the operator to get the starting value
    expression->operator()(this, field);
    m_assert(field->ghost_status(), "the ghost status should be valid here...");

    // refine given the value
    Field*    myfield = field();
    SetValue* my_expr = expression();
    AdaptMagic(field, nullptr, &cback_StatusCheck, &cback_StatusCheck, reinterpret_cast<void*>(myfield), &cback_ValueFill, reinterpret_cast<void*>(my_expr));

    // we modified one block after another, so we set the ghost value from the SetValue
    field->ghost_status(expression->do_ghost());
    //-------------------------------------------------------------------------
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
void Grid::Adapt(m_ptr<list<Patch> > patches) {
    m_begin;
    //-------------------------------------------------------------------------
    if (patches->size() == 0) {
        return;
    }
    // set the recursive mode to true
    // SetRecursiveAdapt(recursive_adapt);
    // go to the magical adapt function
    // Adapt(nullptr, &cback_Patch, &cback_Patch, reinterpret_cast<void*>(patches()), &cback_AllocateOnly, nullptr);
    AdaptMagic(nullptr, patches, &cback_StatusCheck, &cback_StatusCheck, reinterpret_cast<void*>(patches()), &cback_UpdateDependency, nullptr);
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief Adapt the grid, this function performs the actual adaptation for all the others
 * 
 * @warning This function is a "do all the magic function" and therefore should be changed carefully.
 * 
 * @param field_detail if not empty, the field used to computer the details on the blocks
 * @param patches if not empty, a list of patches to dictate the grid layout
 * @param coarsen_cback p4est callback function for deciding on coarsening
 * @param refine_cback p4est callback function for deciding on refinement
 * @param coarseref_cback_ptr pointer to data used in these refine/coarsen p4est callback functions
 * @param interpolate_fct p4est callback function for interpolation: do the actual refinement/coarsening
 * @param interpolate_ptr pointer to data used in the interpolate callback function
 */
void Grid::AdaptMagic(/* criterion */ m_ptr<Field> field_detail, m_ptr<list<Patch> > patches,
                      /* p4est coarsen/refine */ cback_coarsen_citerion_t coarsen_cback, cback_refine_criterion_t refine_cback, void* coarseref_cback_ptr,
                      /* p4est interpolate */ cback_interpolate_t interpolate_fct, void* interpolate_ptr) {
    m_begin;
    m_assert(interpolate_fct != nullptr, "the interpolation function cannot be a nullptr");
    m_assert(cback_criterion_ptr_ == nullptr, "the pointer `cback_criterion_ptr` must be  null");
    m_assert(cback_interpolate_ptr_ == nullptr, "the pointer `cback_interpolate_ptr` must be  null");
    m_assert(p4est_forest_->user_pointer == nullptr, "we must reset the user_pointer to null");
    m_assert(field_detail.IsEmpty() || patches.IsEmpty(), "you cannot give both a field for detail computation and a patch list");
    // m_assert(!(recursive_adapt() && interp_fct == &cback_UpdateDependency), "we cannot use the update dependency in a recursive mode");
    //-------------------------------------------------------------------------
    m_profStart(prof_, "adaptation");
    //................................................
    // pre-assigne the profiling in case every cpu doesn't enter it
    m_profInitLeave(prof_, "criterion detail");
    m_profInitLeave(prof_, "solve dependency");

    // inform we start the mess
    m_log("--> grid adaptation started... (recursive = %d)", recursive_adapt());

    //................................................
    DoOpMesh(nullptr, &GridBlock::ResetStatus, this);

    //................................................
    // store the ptrs and the grid
    cback_criterion_ptr_        = coarseref_cback_ptr;
    cback_interpolate_ptr_      = interpolate_ptr;
    p4est_forest_->user_pointer = reinterpret_cast<void*>(this);

    //................................................
    // coarsening, need to have the 8 children on the same rank
    iter_t   iteration              = 0;
    iblock_t global_n_quad_to_adapt = 0;
    do {
        // reset the adapt counter
        n_quad_to_adapt_       = 0;
        global_n_quad_to_adapt = 0;

        //................................................

        // if we consider a field, compute the criterion based on the field
        // if not, just go ahead
        if (!field_detail.IsEmpty()) {
            // compute the details
            m_log("using details");
            DoOpTree(nullptr, &GridBlock::UpdateStatusFromCriterion, this,
                     m_ptr<const Wavelet>(interp_), rtol_, ctol_, m_ptr<const Field>(field_detail), m_ptr<Prof>(prof_));
        }
        if (!patches.IsEmpty()) {
            // get the patches processed
            m_log("using patches");
            DoOpTree(nullptr, &GridBlock::UpdateStatusFromPatches, this,
                     m_ptr<const Wavelet>(interp_), patches, m_ptr<Prof>(prof_));
        }

        // // need to put the status in the array
        // DoOpTree(nullptr, &GridBlock::SetNewByCoarsening, this, m_ptr<bool>(coarsen_status_));

        // // need to know now otherwise I loose the access
        // ExchangeStatus_PostStart_();
        // ExchangeStatus_CompleteWait_();

        //................................................
        // after this point, we cannot access the old blocks anymore, p4est will destroy the access.
        // we still save them as dependencies but all the rest is gone.
        //................................................
        // need to update on the neighbors before doing anything
        DestroyMeshGhost();
        DestroyAdapt();

        // coarsening for p4est-> only one level
        // The limit in levels are handled directly on the block, not in p4est
        if (coarsen_cback != nullptr) {
            m_profStart(prof_, "p4est coarsen");
            p8est_coarsen_ext(p4est_forest_, false, 0, coarsen_cback, nullptr, interpolate_fct);
            m_profStop(prof_, "p4est coarsen");
        }

        // refinement -> only one level
        // The limit in levels are handled directly on the block, not in p4est
        if (refine_cback != nullptr) {
            m_profStart(prof_, "p4est refine");
            p8est_refine_ext(p4est_forest_, false, P8EST_QMAXLEVEL, refine_cback, nullptr, interpolate_fct);
            m_profStop(prof_, "p4est refine");
        }

        // if we are recursive, we need to update the rank partitioning and check for new block to change
        if (recursive_adapt()) {
            // if we balance the grid, it enters a endless circle as the coarsening and refinement corrects the balanced parition
            // no balancing on the grid but corrects the proc distribution to have everybody on the same cpu
            m_profStart(prof_, "partition init");
            Partitioner partition(&fields_, this, true);
            m_profStop(prof_, "partition init");
            m_profStart(prof_, "partition comm");
            partition.Start(&fields_, M_FORWARD);
            partition.End(&fields_, M_FORWARD);
            m_profStop(prof_, "partition comm");
        }

        // sum over the ranks and see if we keep going
        m_assert(n_quad_to_adapt_ < std::numeric_limits<int>::max(), "we must be smaller than the integer limit");
        MPI_Allreduce(&n_quad_to_adapt_, &global_n_quad_to_adapt, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        m_log("we have %d blocks to adapt", global_n_quad_to_adapt);

        // increment the counter
        ++iteration;
    } while (global_n_quad_to_adapt > 0 && recursive_adapt() && iteration < P8EST_QMAXLEVEL);

    //................................................
    // get the 2:1 constrain on the grid
    m_profStart(prof_, "p4est balance");
    p8est_balance_ext(p4est_forest_, P8EST_CONNECT_FULL, nullptr, interpolate_fct);
    m_profStop(prof_, "p4est balance");

    //................................................
    // exchange the status
    // ExchangeStatus_PostStart_();

    // Solve the dependencies if some have been created -> no dep created if we are recursive
    // this is check in the callback UpdateDependency function
    if (!recursive_adapt()) {
        // warn the user that we do not interpolate a temporary field
        for (auto fid = FieldBegin(); fid != FieldEnd(); ++fid) {
            m_log("field <%s> %s", fid->second->name().c_str(), fid->second->is_temp() ? "is discarded" : "will be interpolated");
        }
        // do the solve dependency, level by level!
        m_verb("solve dependencies");
        DoOpTree(nullptr, &GridBlock::SolveDependency, this, interp_, FieldBegin(), FieldEnd(), prof_);

        // we need to process
    } else {
        m_log("no depency solving as we were recursivelly refining");
    }

    //................................................
    // finally fix the rank partition
    m_profStart(prof_, "partition init");
    Partitioner partition(&fields_, this, true);
    m_profStop(prof_, "partition init");
    m_profStart(prof_, "partition comm");
    partition.Start(&fields_, M_FORWARD);
    partition.End(&fields_, M_FORWARD);
    m_profStop(prof_, "partition comm");

    //................................................
    // create a new ghost and mesh as the partioning is done
    SetupMeshGhost();
    SetupAdapt();

    //................................................
    ghost_->UpdateStatus();

    // solve resolution jump if needed
    m_verb("solve jump resolution");
    DoOpTree(nullptr, &GridBlock::SolveResolutionJump, this, interp_, FieldBegin(), FieldEnd(), prof_);

    //................................................
    // set the ghosting fields as non-valid
    for (auto fid = fields_.begin(); fid != fields_.end(); fid++) {
        fid->second->ghost_status(false);
    }

    //................................................
    // reset the forest pointer
    cback_criterion_ptr_        = nullptr;
    cback_interpolate_ptr_      = nullptr;
    p4est_forest_->user_pointer = nullptr;

    // reset the recursive to false
    SetRecursiveAdapt(false);

    // and finally end that stuff
    m_profStop(prof_, "adaptation");
    //-------------------------------------------------------------------------
    m_assert(cback_criterion_ptr_ == nullptr, "the pointer `cback_criterion_ptr` must be  null");
    m_assert(cback_interpolate_ptr_ == nullptr, "the pointer `cback_interpolate_ptr` must be  null");
    m_assert(p4est_forest_->user_pointer == nullptr, "we must reset the user_pointer to null");
    level_t min_level = this->MinLevel();
    level_t max_level = this->MaxLevel();
    m_log("--> grid adaptation done: now %ld blocks on %ld trees using %d ranks and %d threads (level from %d to %d)", p4est_forest_->global_num_quadrants, p4est_forest_->trees->elem_count, p4est_forest_->mpisize, omp_get_max_threads(), min_level, max_level);
    m_end;
}

void Grid::DumpDetails(m_ptr<Field> criterion, m_ptr<Field> details) {
    m_begin;
    m_assert(criterion->lda() == details->lda(), "field <%s> and <%s> must have the same size", criterion->name().c_str(), details->name().c_str());
    //-------------------------------------------------------------------------
    this->GhostPull(criterion);

    // const Wavelet* const_interp    = interp_;
    // const Field*   const_criterion = criterion;
    // const Field*   const_details   = details;
    DoOpMesh(nullptr, &GridBlock::ComputeDetails, this, interp(), criterion, details);
    details->ghost_status(false);
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief Compute the wavelet transform of a current block and it's refinement/coarsening status
 * 
 * a block is refined if needed
 * a block is coarsened if possible and if none of its finer neighbors are coarsened
 * 
 * @param field 
 */
void Grid::ExchangeStatus_PostStart_() const {
    // m_assert(coarsen_status_ != nullptr, "we need the ghosting to be alive");
    // m_assert(coarsen_status_window_ != MPI_WIN_NULL, "the window cannot be null");
    // m_assert(ghost_ != nullptr, "the ghost must exist");
    // m_assert(ghost_->mirror_origin_group() != MPI_GROUP_NULL, "call the InitComm function first!");
    // m_assert(ghost_->mirror_target_group() != MPI_GROUP_NULL, "call the InitComm function first!");
    // // m_assert(field->ghost_status(), "the field <%s> must have up to date ghost values", field->name().c_str());
    // //-------------------------------------------------------------------------

   
    //-------------------------------------------------------------------------
}
void Grid::ExchangeStatus_CompleteWait_() const {
    // m_assert(coarsen_status_ != nullptr, "we need the ghosting to be alive");
    // m_assert(coarsen_status_window_ != MPI_WIN_NULL, "the window cannot be null");
    // m_assert(ghost_ != nullptr, "the ghost must exist");
    // m_assert(ghost_ != nullptr, "we need the ghosting to be alive");
    // m_assert(field->ghost_status(), "the field <%s> must have up to date ghost values", field->name().c_str());
    //-------------------------------------------------------------------------
    // // close the access epochs
    // if (ghost_->mirror_target_group() != MPI_GROUP_EMPTY) {
    //     MPI_Win_complete(coarsen_status_window_);
    // }
    // // close the exposure epochs
    // MPI_Win_wait(coarsen_status_window_);

    //-------------------------------------------------------------------------
}