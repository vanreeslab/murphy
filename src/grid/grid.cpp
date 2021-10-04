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
Grid::Grid(const level_t ilvl, const bool isper[3], const lid_t l[3], BlockDataType block_type, MPI_Comm comm, Prof* const prof)
    : ForestGrid(ilvl, isper, l, block_type, BlockPointerSize(block_type), comm) {
    m_begin;
    //-------------------------------------------------------------------------
    // profiler
    prof_ = prof;
    // create a default Wavelet -> default is M_WAVELET_N and M_WAVELET_NT
    interp_ = new InterpolatingWavelet();

    // create the associated blocks
    p8est_iter_volume_t callback_create = get_cback_CreateBlock(block_type);
    p8est_iterate(p4est_forest_, nullptr, nullptr, callback_create, nullptr, nullptr, nullptr);

    // partition the grid to have compatible grid
    Partitioner part = Partitioner(&fields_, this, true);
    // part.Start(&fields_, M_FORWARD);
    // part.End(&fields_, M_FORWARD);
    part.SendRecv(&fields_, M_FORWARD);

    // setup the ghost stuctures as the mesh will not change anymore
    SetupMeshGhost();
    //-------------------------------------------------------------------------
    m_log("uniform grid created with %ld blocks on %ld trees using %d ranks and %d threads", p4est_forest_->global_num_quadrants, p4est_forest_->trees->elem_count, p4est_forest_->mpisize, omp_get_max_threads());
    m_end;
}

// /**
//  * @brief Copy the ForestGrid part from a grid
//  *
//  * @param grid the source grid
//  */
// void Grid::CopyFrom(const Grid* grid) {
//     m_begin;
//     //-------------------------------------------------------------------------
//     this->ForestGrid::CopyFrom(grid);
//     // copy the field mapping
//     for (auto iter = grid->FieldBegin(); iter != grid->FieldEnd(); iter++) {
//         string name   = iter->first;
//         Field* fid    = iter->second;
//         fields_[name] = fid;
//     }
//     // copy the profiler
//     prof_ = grid->profiler();
//     //-------------------------------------------------------------------------
//     m_end;
// }

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
    // destroy the remaining blocks
    if (is_connect_owned_) {
        p8est_iterate(p4est_forest_, nullptr, nullptr, cback_DestroyBlock, nullptr, nullptr, nullptr);
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
    m_assert(this->MinLevel() >= min, "trying to impose a min level = %d while blocks exist on a lower one = %d", min, this->MinLevel());
    m_assert(this->MaxLevel() <= max, "trying to impose a max level = %d while blocks exist on a higher one = %d", max, this->MaxLevel());
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
        delete ghost_;
        ghost_ = nullptr;
    }
    this->DestroyP4estMeshAndGhost();
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
    for (const auto fid : fields_) {
        memsize += p4est_forest_->local_num_quadrants * (M_N * M_N * M_N) * fid.second->lda() * sizeof(real_t);
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

bool Grid::IsAField(const Field* const field) const {
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
void Grid::AddField(Field* field) {
    m_begin;
    // m_assert(!field.IsOwned(), "The field cannot be owned as it has not been created here");
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
void Grid::DeleteField(const Field* field) {
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
void Grid::ResetFields(const std::map<string, Field*>* fields) {
    m_begin;
    //-------------------------------------------------------------------------
    // clear the current map
    fields_.clear();
    // copy the new one
    for (auto fid : *fields) {
        fields_[fid.first] = fid.second;

        // check if we satisfy the requirements on the key
        m_assert(fid.first == fid.second->name(), "the key of the map must be the name of the field");
    }
    //-------------------------------------------------------------------------
    m_end;
}

void Grid::GhostPull_SetLength(const Field* field, bidx_t ghost_len[2]) const {
    m_begin;
    m_assert(!(field == nullptr), "the source field cannot be empty");
    m_assert(IsAField(field), "the field does not belong to this grid");
    m_assert(ghost_ != nullptr, "The ghost structure is not valid, unable to use it");
    //-------------------------------------------------------------------------
    m_log("ghost check: field <%s> is %s (requested %d %d, provided %d %d)", field->name().c_str(), field->ghost_status(ghost_len) ? "OK" : "to be computed", ghost_len[0], ghost_len[1], field->get_ghost_len(0), field->get_ghost_len(1));
    if (!field->ghost_status(ghost_len)) {
        ghost_->SetLength(ghost_len);
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
void Grid::GhostPull_Post(const Field* field, const sid_t ida, const bidx_t ghost_len[2]) const {
    m_begin;
    m_assert(0 <= ida && ida < field->lda(), "the ida is not within the field's limit");
    m_assert(!(field == nullptr), "the source field cannot be empty");
    m_assert(IsAField(field), "the field does not belong to this grid");
    m_assert(ghost_ != nullptr, "The ghost structure is not valid, unable to use it");
    //-------------------------------------------------------------------------
    if (!field->ghost_status(ghost_len)) {
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
void Grid::GhostPull_Wait(const Field* field, const sid_t ida, const bidx_t ghost_len[2]) const {
    m_begin;
    m_assert(0 <= ida && ida < field->lda(), "the ida is not within the field's limit");
    m_assert(!(field == nullptr), "the source field cannot be empty");
    m_assert(IsAField(field), "the field does not belong to this grid");
    m_assert(ghost_ != nullptr, "The ghost structure is not valid, unable to use it");
    //-------------------------------------------------------------------------
    if (!field->ghost_status(ghost_len)) {
        ghost_->PullGhost_Wait(field, ida);
    }
    //-------------------------------------------------------------------------
    m_end;
}

void Grid::GhostPull(Field* field, const BlockOperator* op) const {
    m_begin;
    //-------------------------------------------------------------------------
    bidx_t ghost_len[2];
    op->GhostLengthNeed(ghost_len);
    GhostPull(field, ghost_len);
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief Pull the ghost points (take the values from the neighbors): after this function, the ghost status of the field is set as up-to-date
 * 
 * @param field the field which requires the ghost
 */
void Grid::GhostPull(Field* field, const bidx_t ghost_len_usr[2]) const {
    m_begin;
    m_assert(!(field == nullptr), "the source field cannot be empty");
    m_assert(ghost_ != nullptr, "The ghost structure is not valid, unable to use it");
    //-------------------------------------------------------------------------
    // get the real ghost length
    bidx_t ghost_len[2] = {ghost_len_usr[0], ghost_len_usr[1]};
    GhostPull_SetLength(field, ghost_len);

    // start the send in the first dimension
    // m_log("ghost check: field <%s> is %s", field->name().c_str(), field->ghost_status(ghost_len) ? "OK" : "to be computed");
    m_profStart(prof_, "pull ghost");
    for (lda_t ida = 0; ida < field->lda(); ++ida) {
        m_verb("ghosting post field <%s> in dir %d", field->name().c_str(), ida);
        GhostPull_Post(field, ida, ghost_len);
        m_verb("ghosting wait field <%s> in dir %d", field->name().c_str(), ida);
        GhostPull_Wait(field, ida, ghost_len);
    }
    m_profStop(prof_, "pull ghost");
    // set that everything is ready for the field
    const bidx_t ghost_len_actual[2] = {m_max(ghost_len[0], field->get_ghost_len(0)),
                                        m_max(ghost_len[1], field->get_ghost_len(1))};
    field->ghost_len(ghost_len_actual);
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
    m_assert(refine_tol > 0.0, "The refinement tolerance = %e must be > 0.0", refine_tol);
    m_assert(coarsen_tol >= 0.0, "The coarsening tolerance = %e must be >= 0.0", coarsen_tol);
    //-------------------------------------------------------------------------
    // if ctol is smaller than epsilon, just remember epsilon
    ctol_ = m_max(coarsen_tol, std::numeric_limits<real_t>::epsilon());
    rtol_ = refine_tol;
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
void Grid::Refine(Field* field) {
    m_begin;
    m_assert(IsAField(field), "the field must already exist on the grid!");
    m_assert(!recursive_adapt(), "we cannot refine recursivelly here");
    //-------------------------------------------------------------------------
    // get the ghost length
    const bidx_t ghost_len[2] = {interp_->nghost_front(), interp_->nghost_back()};
    // compute the ghost needed by the interpolation of every other field in the grid
    for (auto fid : fields_) {
        Field* cur_field = fid.second;
        if (!cur_field->is_temp()) {
            GhostPull(cur_field, ghost_len);
        }
    }

    const cback_interpolate_t cback_update_dependency = get_cback_UpdateDependency(block_type_);
    AdaptMagic(field, nullptr, nullptr, &cback_StatusCheck, nullptr, cback_update_dependency, nullptr);
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
void Grid::Coarsen(Field* field) {
    m_begin;
    m_assert(IsAField(field), "the field must already exist on the grid!");
    // m_assert(!recursive_adapt(), "we cannot refine recursivelly here");
    //-------------------------------------------------------------------------
    // // get the ghost length
    // const bidx_t ghost_len[2] = {interp_->nghost_front(), interp_->nghost_back()};
    // // compute the ghost needed by the interpolation of every other field in the grid
    // for (auto* fid : fields_) {
    //     Field* cur_field = fid.second;
    //     if (!cur_field->is_temp()) {
    //         GhostPull(cur_field, ghost_len);
    //     }
    // }

    const cback_interpolate_t cback_update_dependency = get_cback_UpdateDependency(block_type_);
    AdaptMagic(field, nullptr, &cback_StatusCheck, nullptr, nullptr, cback_update_dependency, nullptr);
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief Adapt the grid based on field which is used as a criterion
 * 
 * The adaptation might be recursive, see Grid::SetRecursiveAdapt() function
 * 
 * @param field the field used for the criterion
 */
void Grid::Adapt(Field* field) {
    m_begin;
    m_assert(IsAField(field), "the field must already exist on the grid!");
    //-------------------------------------------------------------------------
    const cback_interpolate_t cback_update_dependency = get_cback_UpdateDependency(block_type_);
    AdaptMagic(field, nullptr, &cback_StatusCheck, &cback_StatusCheck, nullptr, cback_update_dependency, nullptr);
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief Adapt the grid given an analytical expression for the designated field (typically an initial condition)
 * 
 * The anaylical expression is used to start and then is used everytime a new block is created.
 * Other fields than field might be garbage
 * 
 * Notes: 
 * - The adaptation can be recursive
 * - if the ghost provided by the SetValue are not sufficient, the ghosts will be recomputed
 * 
 * @param field the field to adapt
 * @param expr the analytical expression to use as a SetValue object
 */
void Grid::Adapt(Field* field, const SetValue* expr) {
    m_begin;
    m_assert(IsAField(field), "the field must already exist on the grid!");
    //-------------------------------------------------------------------------
    // apply the operator to get the starting value
    (*expr)(this, field);

    // refine given the value, have to remove the cast to allow the cas in void* (only possible way, sorry)
    void* my_expr = static_cast<void*>(const_cast<SetValue*>(expr));
    const cback_interpolate_t cback_value_fill = get_cback_ValueFill(block_type_);
    AdaptMagic(field, nullptr, &cback_StatusCheck, &cback_StatusCheck, static_cast<void*>(field), cback_value_fill, my_expr);

    // apply the operator to get the starting value
    // it gets rid of whatever smoothing has been done
    // (*expr)(this, field);
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
void Grid::Adapt(list<Patch>* patches) {
    m_begin;
    //-------------------------------------------------------------------------
    if (patches->empty()) {
        return;
    }
    const cback_interpolate_t cback_update_dependency = get_cback_UpdateDependency(block_type_);
    AdaptMagic(nullptr, patches, &cback_StatusCheck, &cback_StatusCheck, nullptr, cback_update_dependency, nullptr);
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
 * @param max_detail stores the value of the max_detail per dimension of the field, if not nullptr
 */
void Grid::AdaptMagic(/* criterion */ Field* field_detail, list<Patch>* patches,
                      /* p4est coarsen/refine */ cback_coarsen_citerion_t coarsen_cback, cback_refine_criterion_t refine_cback, void* coarseref_cback_ptr,
                      /* p4est interpolate */ cback_interpolate_t interpolate_fct, void* interpolate_ptr) {
    m_begin;
    m_assert(interpolate_fct != nullptr, "the interpolation function cannot be a nullptr");
    m_assert(cback_criterion_ptr_ == nullptr, "the pointer `cback_criterion_ptr` must be  nullptr");
    m_assert(cback_interpolate_ptr_ == nullptr, "the pointer `cback_interpolate_ptr` must be  nullptr");
    m_assert(p4est_forest_->user_pointer == nullptr, "we must reset the user_pointer to nullptr");
    m_assert((field_detail == nullptr) || (patches == nullptr), "you cannot give both a field for detail computation and a patch list");
    // m_assert(!(recursive_adapt() && interp_fct == get_cback_UpdateDependency(block_type_), "we cannot use the update dependency in a recursive mode");
    //-------------------------------------------------------------------------
    m_profStart(prof_, "adaptation");
    //................................................
    // pre-assigne the profiling in case every cpu doesn't enter it
    // m_profInitLeave(prof_, "criterion detail");
    // m_profInitLeave(prof_, "patch");
    // m_profInitLeave(prof_, "solve dependency");
    // m_profInitLeave(prof_, "smooth jump");

    // inform we start the mess
    string msg = "grid adaptation started... (recursive = %d, ";
    if (!(field_detail == nullptr)) {
        msg += "using details";
    }
    if (!(patches == nullptr)) {
        msg += "using patches";
    }
    msg += ") -> %d fields";
    m_log(msg.c_str(), recursive_adapt(), fields_.size());

    //................................................
    // store the ptrs and the grid
    cback_criterion_ptr_        = coarseref_cback_ptr;
    cback_interpolate_ptr_      = interpolate_ptr;
    p4est_forest_->user_pointer = reinterpret_cast<void*>(this);

    // reset the status of everyblock to a neutral one
    m_profStart(prof_, "reset");
    DoOpMesh(nullptr, &GridBlock::StatusReset, this);
    m_profStop(prof_, "reset");

    //................................................
    // coarsening, need to have the 8 children on the same rank
    iter_t   iteration              = 0;
    iblock_t global_n_quad_to_adapt = 0;

    do {
        m_log("----> adaptation iteration %d", iteration);
        m_log_level_plus;
        //......................................................................
        // cleanup the status
        DoOpMesh(nullptr, &GridBlock::StatusCleanup, this);
#ifndef NDEBUG
        {
            const level_t min_level = this->MinLevel();
            const level_t max_level = this->MaxLevel();
            m_log("--> iteration on a grid with %ld blocks on %ld trees using %d ranks and %d threads (level from %d to %d)", p4est_forest_->global_num_quadrants, p4est_forest_->trees->elem_count, p4est_forest_->mpisize, omp_get_max_threads(), min_level, max_level);
        }
#endif

        //......................................................................
        // get the desired status
        if (!(field_detail == nullptr)) {
            m_assert(!field_detail->is_temp(), "The criterion field cannot be temporary");
            // ghost the dimension 0
            m_profStart(prof_, "ghost for criterion");
            
            // get the length
            bidx_t ghost_len[2] = {interp_->nghost_front(), interp_->nghost_back()};
            GhostPull_SetLength(field_detail, ghost_len);

            // start the ghosting in the first dimension only
            GhostPull_Post(field_detail, 0, ghost_len);
            GhostPull_Wait(field_detail, 0, ghost_len);
            m_profStop(prof_, "ghost for criterion");
            // now we can start the next dimensions
            for (lda_t ida = 1; ida < field_detail->lda(); ++ida) {
                // start the ghosts for the next dimension
                m_profStart(prof_, "ghost for criterion");
                GhostPull_Post(field_detail, ida, ghost_len);
                m_profStop(prof_, "ghost for criterion");

                // compute the criterion on the previous dimension
                m_profStart(prof_, "criterion");
                DoOpMesh(nullptr, &GridBlock::UpdateStatusFromCriterion, this, interp_, rtol_, ctol_, field_detail, ida - 1);
                m_profStop(prof_, "criterion");

                // finish the ghost for the current dimension
                m_profStart(prof_, "ghost for criterion");
                GhostPull_Wait(field_detail, ida, ghost_len);
                m_profStop(prof_, "ghost for criterion");
            }
            // finally, compute on the last dimension of the field
            m_profStart(prof_, "criterion");
            const lda_t criterion_dim = field_detail->lda() - 1;
            DoOpMesh(nullptr, &GridBlock::UpdateStatusFromCriterion, this, interp_, rtol_, ctol_, field_detail, criterion_dim);
            m_profStop(prof_, "criterion");

            // register the computed ghosts
            field_detail->ghost_len(ghost_len);
        }
        // if needed, compute some patches
        if (!(patches == nullptr)) {
            // get the patches processed
            m_profStart(prof_, "patch");
            DoOpMesh(nullptr, &GridBlock::UpdateStatusFromPatches, this, interp_, patches);
            m_profStop(prof_, "patch");
        }
        //......................................................................
#ifndef NDEBUG
        if (field_detail != nullptr) {
            // check/get the max detail on the current grid
            real_t det_maxmin[2];
            this->MaxMinDetails(field_detail, det_maxmin);
            m_log("--> before adaptation: rtol = %e, max detail = %e", this->rtol(), det_maxmin[0]);
        }
#endif

        //......................................................................
        // every block has its own desired status by now: finer, coarser or same
        // we need to apply a few policies to make sure that
        // (1) two neighboring blocks don't do the opposite operation (coarsen and refine)
        // (2)
        {
            m_profStart(prof_, "update status");
            // first apply the local policies (max/min levels)
            DoOpMesh(nullptr, &GridBlock::UpdateStatusFromLevel, this, level_limit_min_, level_limit_max_);

            // allocate the neighbor status vectors on everyblock
            ghost_->SyncStatusInit();

            // I need to make sure that the if one block needs to refine, all the coarser neighbors get that signal
            // first step is to put refine to everybody
            const level_t max_level = MaxLevel();
            const level_t min_level = MinLevel();
            for (level_t il = (max_level - 1); il >= min_level; --il) {
                // the finer level fills the status
                ghost_->SyncStatusFill(il + 1, il + 1);
                // get the update in the current level
                ghost_->SyncStatusUpdate(il, il);
                // forward the refinement criterion: if my finer neighbor has to refine, I refine as well
                DoOpMeshLevel(nullptr, &GridBlock::UpdateStatusForwardRefinement, this, il);
            }

            // now adapt my status based on the status of my neighbors
            ghost_->SyncStatusFill();
            ghost_->SyncStatusUpdate();
            DoOpMesh(nullptr, &GridBlock::UpdateStatusFromGlobalPolicy, this);

            ghost_->SyncStatusFinalize();
            m_profStop(prof_, "update status");
        }

        //................................................
        // after this point, we cannot access the old blocks anymore, p4est will destroy the access.
        // we still save them as dependencies but all the rest is gone.
        m_profStart(prof_, "destroy mesh and ghost");
        DestroyMeshGhost();
        m_profStop(prof_, "destroy mesh and ghost");

        //................................................
        // reset the adapt counter, it will be updated in the interpolate fct
        n_quad_to_refine_  = 0;
        n_quad_to_coarsen_ = 0;
        // WARNING: always try to COARSEN first (no idea why but the other way around doesn't work!)
        // coarsening for p4est-> only one level
        // The limit in levels are handled directly on the block, not in p4est
        if (coarsen_cback != nullptr) {
            m_profStart(prof_, "p4est coarsen");
            p8est_coarsen_ext(p4est_forest_, 0, 0, coarsen_cback, nullptr, interpolate_fct);
            m_profStop(prof_, "p4est coarsen");
        }
        // refinement -> only one level
        // The limit in levels are handled directly on the block, not in p4est
        if (refine_cback != nullptr) {
            m_profStart(prof_, "p4est refine");
            p8est_refine_ext(p4est_forest_, 0, P8EST_QMAXLEVEL, refine_cback, nullptr, interpolate_fct);
            m_profStop(prof_, "p4est refine");
        }

        // get the 2:1 constrain on the grid, should be guaranteed by the criterion, but just in case
        // m_profStart(prof_, "p4est balance");
        // p8est_balance_ext(p4est_forest_, P8EST_CONNECT_FULL, nullptr, interpolate_fct);
        // m_profStop(prof_, "p4est balance");

        //................................................
        // solve the dependencies on the grid
        // warn the user that we do not interpolate a temporary field
        for (auto fid = FieldBegin(); fid != FieldEnd(); ++fid) {
            m_log("field <%s> %s", fid->second->name().c_str(), fid->second->is_temp() ? "is discarded" : "will be interpolated");
        }
        m_profStart(prof_, "solve dependency");
        DoOpTree(nullptr, &GridBlock::SolveDependency, this, interp_, FieldBegin(), FieldEnd());
        m_profStop(prof_, "solve dependency");

        //................................................
        // update the rank partitioning and check for new block to change
        m_profStart(prof_, "partition init");
        Partitioner partition(&fields_, this, true);
        m_profStop(prof_, "partition init");
        m_profStart(prof_, "partition comm");
        partition.SendRecv(&fields_, M_FORWARD);
        // partition.Start(&fields_, M_FORWARD);
        // partition.End(&fields_, M_FORWARD);
        m_profStop(prof_, "partition comm");

        //................................................
        // create a new ghost and mesh as the partioning is done
        m_profStart(prof_, "setup mesh and ghost");
        SetupMeshGhost();
        m_profStop(prof_, "setup mesh and ghost");

        //................................................
        // solve the jump in resolution
        m_profStart(prof_, "update status");
        // get to know the status of my neighbors
        ghost_->SyncStatusInit();
        ghost_->SyncStatusFill();
        ghost_->SyncStatusUpdate();
        // ghost_->SyncStatusFinalize();
        m_profStop(prof_, "update status");

        // solve resolution jump if needed
        m_verb("solve jump resolution");
        m_profStart(prof_, "smooth jump");
        DoOpTree(nullptr, &GridBlock::SmoothResolutionJump, this, interp_, FieldBegin(), FieldEnd());
        m_profStop(prof_, "smooth jump");

        ghost_->SyncStatusFinalize();

        //................................................
        // sum over the ranks and see if we keep going
        m_assert(n_quad_to_coarsen_ < std::numeric_limits<int>::max(), "we must be smaller than the integer limit");
        m_assert(n_quad_to_refine_ < std::numeric_limits<int>::max(), "we must be smaller than the integer limit");
        m_assert((n_quad_to_coarsen_ % 8) == 0, "the number of quad to coarsen must be a multiple of 8 isntead of %d", n_quad_to_coarsen_);

        int global_n_adapt[2];
        int n_adapt[2] = {n_quad_to_coarsen_, n_quad_to_refine_};
        MPI_Allreduce(n_adapt, global_n_adapt, 2, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        global_n_quad_to_adapt = global_n_adapt[0] + global_n_adapt[1];
        m_log("we have coarsened %d blocks and refined %d blocks -> %d new blocks", global_n_adapt[0], global_n_adapt[1], global_n_adapt[0] / 8 + 8 * global_n_adapt[1]);

        // if we adapted some blocks, then the ghosting is not valid
        if (global_n_quad_to_adapt > 0) {
            for (auto fid : fields_) {
                m_log("changing ghost status of <%s>", fid.second->name().c_str());
                const bidx_t ghost_len[2] = {0, 0};
                fid.second->ghost_len(ghost_len);
            }
        }

#ifndef NDEBUG
        if (field_detail != nullptr) {
            // check/get the max detail on the current grid
            real_t det_maxmin[2];
            this->MaxMinDetails(field_detail, det_maxmin);
            m_log("--> after adaptation: rtol = %e, max detail = %e", this->rtol(), det_maxmin[0]);
        }
#endif
        m_log_level_minus;
        // increment the counter
        ++iteration;
    } while (global_n_quad_to_adapt > 0 && recursive_adapt() && iteration < P8EST_QMAXLEVEL);

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
    m_assert(cback_criterion_ptr_ == nullptr, "the pointer `cback_criterion_ptr` must be  nullptr");
    m_assert(cback_interpolate_ptr_ == nullptr, "the pointer `cback_interpolate_ptr` must be  nullptr");
    m_assert(p4est_forest_->user_pointer == nullptr, "we must reset the user_pointer to nullptr");
    const level_t min_level   = this->MinLevel();
    const level_t max_level   = this->MaxLevel();
    const MemLayout block_layout(M_LAYOUT_BLOCK, M_GS, M_N);
    const real_t  mem_per_dim = p4est_forest_->global_num_quadrants * block_layout.n_elem * sizeof(real_t) / 1.0e+9;

    // m_log_level_minus;
    m_log("--> grid adaptation done: now %ld blocks (%.2e Gb/dim) on %ld trees using %d ranks and %d threads (level from %d to %d)", p4est_forest_->global_num_quadrants, mem_per_dim, p4est_forest_->trees->elem_count, p4est_forest_->mpisize, omp_get_max_threads(), min_level, max_level);
#ifndef NDEBUG
    if (field_detail != nullptr) {
        // check/get the max detail on the current grid
        real_t det_maxmin[2];
        this->MaxMinDetails(field_detail, det_maxmin);
        m_log("rtol = %e, max detail = %e", this->rtol(), det_maxmin[0]);
        m_assert(this->rtol() >= det_maxmin[0], "the max detail cannot be > than the tol: %e vs %e", det_maxmin[0], this->rtol());
    }
#endif
    m_end;
}

/**
 * @brief Compute and store the details from the criterion field
 * 
 * @param criterion the field on which we compute details
 * @param details the field containing the details
 */
void Grid::StoreDetails(Field* criterion, Field* details) {
    m_begin;
    m_assert(criterion->lda() == details->lda(), "field <%s> and <%s> must have the same size", criterion->name().c_str(), details->name().c_str());
    //-------------------------------------------------------------------------
    {
        const bidx_t ghost_len[2] = {interp_->nghost_front(), interp_->nghost_back()};
        this->GhostPull(criterion, ghost_len);
    }

    DoOpMesh(nullptr, &GridBlock::StoreDetails, this, interp(), criterion, details);
    {
        const bidx_t ghost_len[2] = {0, 0};
        details->ghost_len(ghost_len);
    }
    //-------------------------------------------------------------------------
    m_end;
}

void Grid::MaxMinDetails(Field* criterion, real_t maxmin[2]) {
    //--------------------------------------------------------------------------
    // get the ghosts values
    const bidx_t ghost_len[2] = {interp_->nghost_front(), interp_->nghost_back()};
    this->GhostPull(criterion, ghost_len);

    // get the minmax
    real_t local_maxmin[2];
    local_maxmin[0] = 0.0;
    local_maxmin[1] = std::numeric_limits<real_t>::max();

    DoOpMesh(nullptr, &GridBlock::MaxMinDetails, this, interp(), criterion, local_maxmin, nullptr, 0.0, 0.0, 0);

    // reduce is on the whole mesh
    maxmin[0] = 0.0;
    maxmin[1] = 0.0;
    MPI_Allreduce(local_maxmin, maxmin, 1, M_MPI_REAL, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(local_maxmin + 1, maxmin + 1, 1, M_MPI_REAL, MPI_MIN, MPI_COMM_WORLD);
    //--------------------------------------------------------------------------
}

void Grid::DistributionDetails(const iter_t id, const std::string folder, const std::string suffix, Field* criterion,
                               const short_t n_cat, const real_t max_value_cat) {
    //--------------------------------------------------------------------------
    // get the ghosts values
    const bidx_t ghost_len[2] = {interp_->nghost_front(), interp_->nghost_back()};
    this->GhostPull(criterion, ghost_len);

    bidx_t* n_block_local  = reinterpret_cast<bidx_t*>(m_calloc(n_cat * sizeof(bidx_t)));
    bidx_t* n_block_global = reinterpret_cast<bidx_t*>(m_calloc(n_cat * sizeof(bidx_t)));

    // get the minmax for the whole grid + for each block
    real_t to_trash[2];
    to_trash[0] = 0.0;
    to_trash[1] = std::numeric_limits<real_t>::max();
    DoOpMesh(nullptr, &GridBlock::MaxMinDetails, this, interp(), criterion, to_trash, n_block_local, max_value_cat, 1e-16, n_cat);

    // sum everything on rank 0 to dump
    const rank_t root = 0;
    m_assert(sizeof(bidx_t) == sizeof(int), "the two sizes must match");
    MPI_Reduce(n_block_local, n_block_global, n_cat, MPI_INT, MPI_SUM, root, MPI_COMM_WORLD);

    // dump
    rank_t rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    long nblock = this->global_num_quadrants();

    FILE* file_diag;
    if (rank == root) {
        file_diag = fopen(std::string(folder + "/detail_histogram" + suffix + ".data").c_str(), "a+");
        fprintf(file_diag, "%d;%ld", id, nblock);
        long block_count = 0;
        for (level_t il = 0; il < n_cat; ++il) {
            fprintf(file_diag, ";%d", n_block_global[il]);
            block_count += n_block_global[il];
        }
        fprintf(file_diag, "\n");
        fclose(file_diag);
        m_assert(block_count == nblock, "the two counts must match: %ld %ld", block_count, nblock);
    }

    m_free(n_block_local);
    m_free(n_block_global);

    //--------------------------------------------------------------------------
}