#include "grid.hpp"

#include <omp.h>
#include <p8est_extended.h>

#include "partitioner.hpp"
#include "wavelet/interpolating_wavelet.hpp"

using std::list;
using std::string;
using std::unordered_map;

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
Grid::Grid(const lid_t ilvl, const bool isper[3], const lid_t l[3], MPI_Comm comm, const m_ptr<Prof>& prof)
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
    m_log("ghost check: field <%s> is %s", field->name().c_str(), field->ghost_status() ? "OK" : "to be computed");
    m_profStart(prof_, "pull ghost");
    for (lda_t ida = 0; ida < field->lda(); ++ida) {
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
 * @brief Refine the grid by delta_level locally
 * 
 * @param delta_level the number of level every block will be refined
 */
/**
 * @brief coarsen one field if the criterion matches
 * 
 * @param field the field to coarsen 
 */
void Grid::Refine(m_ptr<Field> field) {
    m_begin;
    m_assert(IsAField(field), "the field must already exist on the grid!");
    m_assert(!recursive_adapt(), "we cannot refine recursivelly here");
    //-------------------------------------------------------------------------
    // compute the ghost needed by the interpolation of every other field in the grid
    for (auto fid = fields_.begin(); fid != fields_.end(); fid++) {
        GhostPull(fid->second);
    }

    // Adapt(reinterpret_cast<void*>(field), nullptr, nullptr, &cback_WaveDetail, &cback_Interpolate);
    // Adapt(reinterpret_cast<void*>(field), nullptr, nullptr, &cback_WaveDetail, &cback_UpdateDependency);
    Adapt(field, nullptr, &cback_StatusCheck, reinterpret_cast<void*>(field()), cback_UpdateDependency, nullptr);

    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief coarsen one field if the criterion matches
 * 
 * @param field the field to coarsen 
 */
void Grid::Coarsen(m_ptr<Field> field) {
    m_begin;
    m_assert(IsAField(field), "the field must already exist on the grid!");
    m_assert(!recursive_adapt(), "we cannot refine recursivelly here");
    //-------------------------------------------------------------------------
    // compute the ghost needed by the interpolation of every other field in the grid
    for (auto fid = fields_.begin(); fid != fields_.end(); fid++) {
        GhostPull(fid->second);
    }

    // Adapt(reinterpret_cast<void*>(field), nullptr, &cback_WaveDetail, nullptr, &cback_Interpolate);
    // Adapt(reinterpret_cast<void*>(field), nullptr, &cback_WaveDetail, nullptr, &cback_UpdateDependency);
    Adapt(field, &cback_StatusCheck, nullptr, reinterpret_cast<void*>(field()), cback_UpdateDependency, nullptr);

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
    Adapt(field, &cback_StatusCheck, &cback_StatusCheck, reinterpret_cast<void*>(field()), cback_UpdateDependency, nullptr);

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
    Adapt(field, &cback_StatusCheck, &cback_StatusCheck, reinterpret_cast<void*>(myfield), &cback_ValueFill, reinterpret_cast<void*>(my_expr));

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
    SetRecursiveAdapt(true);
    // go to the magical adapt function
    Adapt(nullptr, &cback_Patch, &cback_Patch, reinterpret_cast<void*>(patches()), &cback_AllocateOnly, nullptr);
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief 
 * 
 * if the refinement is recursive, I cannot have interp = cback_UpdateDependency
 * 
 * @param field 
 * @param coarsen_crit 
 * @param refine_crit 
 * @param criterion_ptr 
 * @param interp 
 * @param interp_ptr 
 */
void Grid::Adapt(m_ptr<const Field> field, cback_coarsen_citerion_t coarsen_crit, cback_refine_criterion_t refine_crit, void* criterion_ptr, cback_interpolate_t interp, void* interp_ptr) {
    m_begin;
    m_assert(interp != nullptr, "the interpolation function cannot be a nullptr");
    m_assert(cback_criterion_ptr_ == nullptr, "the pointer `cback_criterion_ptr` must be  null");
    m_assert(cback_interpolate_ptr_ == nullptr, "the pointer `cback_interpolate_ptr` must be  null");
    m_assert(p4est_forest_->user_pointer == nullptr, "we must reset the user_pointer to null");
    m_assert(!(recursive_adapt() && interp == &cback_UpdateDependency), "we cannot use the update dependency in a recursive mode");
    //-------------------------------------------------------------------------
    m_profStart(prof_, "adaptation");
    //................................................
    m_profInitLeave(prof_, "criterion detail");
    m_profInitLeave(prof_, "solve dependency");
    // m_profInitLeave(prof_, "p4est coarsen");
    // m_profInitLeave(prof_, "p4est refine");
    // m_profInitLeave(prof_, "p4est balance");
    //................................................
    // log
    m_log("--> grid adaptation started... (recursive = %d)", recursive_adapt());
    //................................................
    // delete the soon-to be outdated ghost and mesh
    DestroyGhost();

    // store the ptrs and the grid
    cback_criterion_ptr_        = criterion_ptr;
    cback_interpolate_ptr_      = interp_ptr;
    p4est_forest_->user_pointer = reinterpret_cast<void*>(this);

    //................................................
    // refinement, can be done recursivelly no matter what the block distributions on the rank
    // the refinement must be done outside otherwise it's an endless circle:
    // refine-coarsen-balance - refine-coarsen-balance - etc

    //................................................
    // coarsening, need to have the 8 children on the same rank
    lid_t iteration              = 0;
    lid_t global_n_quad_to_adapt = 0;
    do {
        // reset the adapt counter
        n_quad_to_adapt_       = 0;
        global_n_quad_to_adapt = 0;

        // if we consider a field, reset the status and comptue the criterion for each block
        if (!field.IsEmpty()) {
            m_log("computing the criterion on field <%s>", field->name().c_str());
            const Wavelet* wavelet_interp = interp_;
            DoOpTree(nullptr, &GridBlock::ResetStatus, this);
            DoOpTree(nullptr, &GridBlock::UpdateStatusCriterion, this, wavelet_interp, rtol_, ctol_, field, prof_);
        }

        // refinement -> only one level
        if (refine_crit != nullptr) {
            m_profStart(prof_, "p4est refine");
            p8est_refine_ext(p4est_forest_, false, P8EST_QMAXLEVEL, refine_crit, nullptr, interp);
            m_profStop(prof_, "p4est refine");
        }

        // coarsening -> only one level
        if (coarsen_crit != nullptr) {
            m_profStart(prof_, "p4est coarsen");
            p8est_coarsen_ext(p4est_forest_, false, 0, coarsen_crit, nullptr, interp);
            m_profStop(prof_, "p4est coarsen");
        }

        // if we are recursive, we need to update the rank partitioning and check for new block to change
        if (recursive_adapt()) {
            // if we balance the grid, it enters a endless circle as the coarsening and refinement corrects the balanced parition
            // no balancing on the grid but corrects the proc distribution
            m_profStart(prof_, "partition init");
            Partitioner partition(&fields_, this, true);
            m_profStop(prof_, "partition init");
            m_profStart(prof_, "partition comm");
            partition.Start(&fields_, M_FORWARD);
            partition.End(&fields_, M_FORWARD);
            m_profStop(prof_, "partition comm");

            // sum over the ranks and see if we keep going
            m_assert(n_quad_to_adapt_ < std::numeric_limits<int>::max(), "we must be smaller than the integer limit");
            MPI_Allreduce(&n_quad_to_adapt_, &global_n_quad_to_adapt, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        }

        m_log("we have %d blocks to adapt locally and %d globally", n_quad_to_adapt_, global_n_quad_to_adapt);

        // increment the counter
        ++iteration;
    } while (global_n_quad_to_adapt != 0 && recursive_adapt() && iteration < P8EST_QMAXLEVEL);

    // get the 2:1 constrain on the grid
    m_profStart(prof_, "p4est balance");
    p8est_balance_ext(p4est_forest_, P8EST_CONNECT_FULL, nullptr, interp);
    m_profStop(prof_, "p4est balance");

    // Solve the dependencies is some have been created -> no dep created if we are recursive
    // this is check in the callback UpdateDependency function
    if (!recursive_adapt()) {
        for (auto fid = FieldBegin(); fid != FieldEnd(); ++fid) {
            m_log("field <%s> %s", fid->second->name().c_str(), fid->second->is_temp() ? "is discarded" : "will be interpolated");
        }
        DoOpTree(nullptr, &GridBlock::SolveDependency, this, interp_, FieldBegin(), FieldEnd(), prof_);
    }

    // finally fix the rank partition
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
    cback_criterion_ptr_        = nullptr;
    cback_interpolate_ptr_      = nullptr;
    p4est_forest_->user_pointer = nullptr;

    // reset the recursive to false
    SetRecursiveAdapt(false);

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