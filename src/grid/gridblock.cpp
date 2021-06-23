#include "gridblock.hpp"

#include <p8est_bits.h>

#include <algorithm>
#include <string>

#include "core/macros.hpp"
#include "p8est_iterate.h"
#include "tools/toolsp4est.hpp"

using std::string;

using GBLocal      = GhostBlock<GridBlock*>;
using GBMirror     = GhostBlock<MPI_Aint>;
using GBPhysic     = PhysBlock;
using ListGBLocal  = std::list<GBLocal*>;
using ListGBMirror = std::list<GBMirror*>;
using listGBPhysic = std::list<GBPhysic*>;

static const lid_t face_start[6][3] = {{0, 0, 0}, {M_N, 0, 0}, {0, 0, 0}, {0, M_N, 0}, {0, 0, 0}, {0, 0, M_N}};

/**
 * @brief return the grid spacing of a coarser block
 * 
 * @param len the physical length of the coarser block
 * @return constexpr real_t the grid spacing in the coarser block
 */
constexpr real_t CoarseHGrid(const real_t len) {
    return len / (real_t)(M_NHALF);
}

/**
 * @brief constructs a new GridBlock given a 3D length and a position
 * 
 * @param length the length of the current block
 * @param xyz the position of the origin, i.e. the left,bottom corner, (x,y,z)
 * @param level the level of the block
 */
GridBlock::GridBlock(const real_t length, const real_t xyz[3], const sid_t level) : CartBlock(length, xyz, level) {
    m_begin;
    //-------------------------------------------------------------------------
    status_lvl_ = M_ADAPT_NONE;

    // init the dependencies
    n_dependency_active_ = 0;
    for (sid_t id = 0; id < P8EST_CHILDREN; ++id) {
        dependency_[id] = nullptr;
    }
    // allocate the coarse ptr
    // size_t alloc_size = m_max(interp->CoarseSize(),);
    // AllocateCoarsePtr(m_blockmemsize(1));
    // m_assert(interp->CoarseSize() )
    coarse_ptr_.Calloc(CartBlockMemNum(1));
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief Destroy the GridBlock
 *
 * delete the ghost ptr if present and delete the associated fields
 *
 */
GridBlock::~GridBlock() {
    //-------------------------------------------------------------------------
    if (coarse_ptr_.IsOwned()) {
        coarse_ptr_.Free();
    }
    // Free the blocks in the mapping
    for (auto it = mem_map_.begin(); it != mem_map_.end(); ++it) {
        it->second.Free();
    }
    //-------------------------------------------------------------------------
}

/**
 * @brief compute the criterion, update the block status accordingly and register the coarsening if needed
 * 
 * The criterion of the block must always be:
 *  rtol >= criterion >= ctol
 * 
 * - if citerion > rtol, we must refine, whatever the other field's dimension value
 * - if citerion < ctol, we can coarsen if this is satisfy by every dimension of the field
 * 
 */
void GridBlock::UpdateStatusFromCriterion(/* params */ const Wavelet*  interp, const real_t rtol, const real_t ctol, const Field*  field_citerion,
                                          /* prof */ Prof*  profiler) {
    //-------------------------------------------------------------------------
    m_assert(rtol > ctol, "the refinement tolerance must be > the coarsening tolerance: %e vs %e", rtol, ctol);
    m_assert(status_lvl_ == M_ADAPT_NONE, "trying to update a status which is already updated");
    m_assert(field_citerion->ghost_status(), "the ghost of <%s> must be up-to-date", field_citerion->name().c_str());
    //-------------------------------------------------------------------------
    m_profStart(profiler, "criterion detail");

    // prevent coarsening if we have finer neighbors
    const bool forbid_coarsening = (local_children_.size() + ghost_children_.size()) > 0;
    const bool forbid_refinement = (local_parent_.size() + ghost_parent_.size()) > 0;

    // I need to visit every dimension and determine if we have to refine and/or coarsen.
    // afterthat we choose given the conservative approach
    bool coarsen = !forbid_coarsening;
    for (lda_t ida = 0; ida < field_citerion->lda(); ida++) {
        // go to the computation
        SubBlock     block_src( -interp->nghost_front(), M_N + interp->nghost_back());
        SubBlock     block_detail( -interp->ndetail_citerion_extend_front(), M_N + interp->ndetail_citerion_extend_back());
        // SubBlock     block_src(this->gs(), this->stride(), -interp->nghost_front(), M_N + interp->nghost_back());
        // SubBlock     block_detail(this->gs(), this->stride(), -interp->ndetail_citerion_extend_front(), M_N + interp->ndetail_citerion_extend_back());

        const real_t norm = interp->Criterion(&block_src, this->data(field_citerion, ida), &block_detail);

        // if the norm is bigger than the refinement tol, we must refine
        bool refine = (norm > rtol) && (!forbid_refinement);
        if (refine) {
            status_lvl_ = M_ADAPT_FINER;
            // finito
            m_profStop(profiler, "criterion detail");
            return;
        }
        // if one dimension is preventing the coarsening, register
        coarsen = coarsen & (norm < ctol);
    }
    // if every field is ok to be coarsened, i.e. the coarsen bool is still true after everything, we coarsen
    // also make sure that we can coarsen!
    status_lvl_ = (coarsen) ? M_ADAPT_COARSER : M_ADAPT_SAME;
    // register the coarsening
    m_profStop(profiler, "criterion detail");
    m_assert(status_lvl_ != M_ADAPT_NONE, "the status of the block cannot be NONE");
    return;
    //-------------------------------------------------------------------------
}

void GridBlock::UpdateStatusFromPatches(/* params */ const Wavelet* interp, std::list<Patch>* patch_list,
                                        /* prof */ Prof* profiler) {
    //-------------------------------------------------------------------------
    m_assert(status_lvl_ == M_ADAPT_NONE, "trying to update a status which is already updated");
    //-------------------------------------------------------------------------

    // prevent coarsening if we have finer neighbors
    const bool forbid_coarsening = (local_children_.size() + ghost_children_.size()) > 0;
    const bool forbid_refinement = (local_parent_.size() + ghost_parent_.size()) > 0;

    // get the block length
    real_t len = p4est_QuadLen(this->level());

    // loop over the patches and determine if I am in it or not
    for (auto patch = patch_list->begin(); patch != patch_list->end(); ++patch) {

        // if we already have the correct level or a higher one, we skip the patch
        if (this->level() > patch->level()) {
            // if not, we have a coarser block and we might want to refine if the location matches
            bool coarsen = !forbid_coarsening;

            for (lda_t id = 0; id < 3; id++) {
                // we have to satisfy both the our max > min and the min < our max
                coarsen = coarsen &&
                          (this->xyz(id) < (patch->origin(id) + patch->length(id))) &&
                          (patch->origin(id) < (this->xyz(id) + len));
            }
            // register the status
            status_lvl_ = (coarsen) ? (M_ADAPT_COARSER) : M_ADAPT_NONE;

            // if we found a matching patch, it's done
            if (coarsen) {
                // m_log("block @ %f %f %f coarsening! ", this->xyz(0), this->xyz(1), this->xyz(2));
                return;
            }

        } else if (this->level() < patch->level()) {
            // if not, we have a coarser block and we might want to refine if the location matches
            bool refine = !forbid_refinement;

            for (lda_t id = 0; id < 3; id++) {
                // we have to satisfy both the our max > min and the min < our max
                refine = refine &&
                         (this->xyz(id) < (patch->origin(id) + patch->length(id))) &&
                         (patch->origin(id) < (this->xyz(id) + len));
            }
            // register the status
            status_lvl_ = (refine) ? (M_ADAPT_FINER) : M_ADAPT_NONE;

            // if we found a matching patch, it's done
            if (refine) {
                // m_log("block @ %f %f %f refining! ", this->xyz(0), this->xyz(1), this->xyz(2));
                return;
            }
        }
    }
    // m_log("block @ %f %f %f not changing! ", this->xyz(0), this->xyz(1), this->xyz(2));
    status_lvl_ = M_ADAPT_SAME;
    m_assert(status_lvl_ != M_ADAPT_NONE, "the status of the block cannot be NONE");
    return;
    //-------------------------------------------------------------------------
}

/**
 * @brief prevents a block from coarsening if one of the neighbor is refining (refine first, always!)
 * 
 */
void GridBlock::UpdateStatusFromPolicy() {
    m_assert(status_lvl_ == M_ADAPT_SAME || status_lvl_ == M_ADAPT_COARSER || status_lvl_ == M_ADAPT_FINER, "the current status must be %d or %d or %d but not %d", M_ADAPT_SAME, M_ADAPT_COARSER, M_ADAPT_FINER, status_lvl_);
    //-------------------------------------------------------------------------
    bool forbid_coarsening = false;

    iblock_t count = 0;
    for (auto gblock : local_parent_) {
        forbid_coarsening = forbid_coarsening || (status_ngh_[count] == M_ADAPT_FINER);
        ++count;
    }
    m_assert(count == local_parent_.size(), "the two numbers must match: %d vs %ld", count, local_parent_.size());
    for (auto gblock : ghost_parent_) {
        forbid_coarsening = forbid_coarsening || (status_ngh_[count] == M_ADAPT_FINER);
        ++count;
    }

    if (forbid_coarsening) {
        StatusAdapt current_status = this->status_level();
        StatusAdapt new_status     = (current_status == M_ADAPT_COARSER) ? M_ADAPT_SAME : current_status;
        this->status_level(new_status);
    }

    m_free(status_ngh_);
    //-------------------------------------------------------------------------
}

void GridBlock::FWTAndGetStatus(const Wavelet*  interp, const real_t rtol, const real_t ctol, const Field*  field_citerion, Prof*  profiler) {
    m_assert(rtol > ctol, "the refinement tolerance must be > the coarsening tolerance: %e vs %e", rtol, ctol);
    m_assert(status_lvl_ == 0, "trying to update a status which is already updated");
    m_assert(field_citerion->ghost_status(), "the ghost of <%s> must be up-to-date", field_citerion->name().c_str());
    //-------------------------------------------------------------------------
    // m_profStart(profiler(), "criterion detail");

    // // I need to visit every dimension and determine if we have to refine and/or coarsen.
    // // afterthat we choose given the conservative approach
    // bool coarsen = true;
    // for (lda_t ida = 0; ida < field_citerion->lda(); ida++) {
    //     // go to the computation
    //     data_ptr     data = this->data(field_citerion, ida);
    //     const real_t norm = interp->FWT(this, data);

    //     // if the norm is bigger than the refinement tol, we must refine
    //     bool refine = norm > rtol;
    //     if (refine) {
    //         status_lvl_ = +1;
    //         // finito
    //         m_profStop(profiler(), "criterion detail");
    //         return;
    //     }
    //     // if one dimension is preventing the coarsening, register
    //     coarsen &= (norm < ctol);
    // }
    // // if every field is ok to be coarsened, i.e. the coarsen bool is still true after everything, we coarsen
    // status_lvl_ = (coarsen) ? -1 : 0;

    // m_profStop(profiler(), "criterion detail");
    m_assert(false, "this function is going to the trash");
    //-------------------------------------------------------------------------
}

/**
 * @brief sets TRUE in the coarsen_vec if the block has been newly created by coarsening
 * 
 * @param qid 
 * @param coarsen_vec 
 */
void GridBlock::SetNewByCoarsening(const qid_t*  qid, short_t* const  coarsen_vec) const {
    // m_assert(status_lvl_ == M_ADAPT_SAME || status_lvl_ == M_ADAPT_NEW_COARSE || status_lvl_ == M_ADAPT_NEW_FINE,"the current status must be %d or %d or %d but not %d", M_ADAPT_SAME, M_ADAPT_NEW_COARSE, M_ADAPT_NEW_FINE, status_lvl_);
    //-------------------------------------------------------------------------
    coarsen_vec[qid->cid] = (short_t)status_lvl_;
    //-------------------------------------------------------------------------
}

/**
 * @brief update my status based on the status from my finer neighbors
 * 
 * allocate the @ref status_siblings_neighbors_ array, which will be destroyed in the @ref SolveNeighbor function.
 */
void GridBlock::GetNewByCoarseningFromNeighbors(const short_t* const  status_vec, MPI_Win status_window) {
    //-------------------------------------------------------------------------
    // get the number of status to obtain
    iblock_t nblocks = (local_parent_.size() + ghost_parent_.size());
    status_ngh_      = reinterpret_cast<short_t*>(m_calloc(nblocks * sizeof(short_t)));

    // loop over the same level neighbors and get the status
    iblock_t count = local_parent_.size();
    for (auto gblock : ghost_parent_) {
        // m_log("count = %d -> requestion block cum id = %ld at rank %d", count, displ, gblock->rank());
        m_assert(sizeof(short_t) == sizeof(short), "the two sizes must match to garantee mpi data types");
        MPI_Get(status_ngh_ + count, 1, MPI_SHORT, gblock->rank(), gblock->cum_block_id(), 1, MPI_SHORT, status_window);
        ++count;
    }
    m_assert(count == nblocks, "the two numbers must match: %d vs %d", count, nblocks);

    // get the local ones now
    count = 0;
    for (auto gblock : local_parent_) {
        // m_log("neighbor id %d gets value %d at id = %d", count, status_vec[gblock->cum_block_id()], gblock->cum_block_id());
        status_ngh_[count] = status_vec[gblock->cum_block_id()];
        ++count;
    }
    m_assert(count == local_parent_.size(), "the two numbers must match: %d vs %ld", count, local_parent_.size());
    //-------------------------------------------------------------------------
}

void GridBlock::SmoothResolutionJump(const Wavelet*  interp, std::map<std::string, Field*  >::const_iterator field_start, std::map<std::string, Field*  >::const_iterator field_end, Prof*  profiler) {
    // the status level has to be 0, otherwise it means that one of the block is not coarsened
    m_assert(status_lvl_ == M_ADAPT_SAME || status_lvl_ == M_ADAPT_NEW_COARSE || status_lvl_ == M_ADAPT_NEW_FINE, "the current status must be %d or %d or %d but not %d", M_ADAPT_SAME, M_ADAPT_NEW_COARSE, M_ADAPT_NEW_FINE, status_lvl_);
    //-------------------------------------------------------------------------
    // reset the temp memory to 0.0
    memset(coarse_ptr_(), 0, CartBlockMemNum(1) * sizeof(real_t));
    data_ptr mask      = coarse_ptr_(0);
    real_t*  mask_data = mask.Write();
    bidx_t   mask_stride = mask.stride();
    
    //................................................
    // lambda to obtain the smoothing pattern
    auto mask_smooth = [=](const iblock_t count, const iface_t ibidule) -> void {
        // create the lambda to put 1.0
        auto set_mask_to_one = [=, &mask_data](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
            mask_data[m_idx(i0, i1, i2, 0, mask_stride)] = 1.0;
        };
        // if the neighbor is a newly created block -> smooth
        if (status_ngh_[count] == M_ADAPT_NEW_COARSE) {
            m_assert(this->status_level() == M_ADAPT_SAME, "if my coarser neighbor has been newly created, I cannot have something different than SAME (now %d)", status_lvl_);

            // get the sign of the ibidule
            real_t sign[3];
            GhostGetSign(ibidule, sign);

            // create the start and end indexes
            bidx_t smooth_start[3], smooth_end[3];
            for (lda_t ida = 0; ida < 3; ++ida) {
                if (sign[ida] > 0.5) {
                    // my ngh assumed 0 details in my block
                    smooth_start[ida] = this->end(ida) - interp->ndetail_citerion_extend_front();
                    // the number of my ngh details influencing my values
                    smooth_end[ida] = this->end(ida) + interp->ndetail_smooth_extend_back();
                } else if (sign[ida] < (-0.5)) {
                    // my ngh assumed 0 details in my block
                    smooth_start[ida] = this->start(ida) - interp->ndetail_smooth_extend_front();
                    // the number of my ngh details influencing my values
                    smooth_end[ida] = this->start(ida) + interp->ndetail_citerion_extend_back();
                } else {
                    // even in the directions orthogonal to ibidule, the details must be killed!
                    // as my neighbor, which might be fine will kill them as well
                    smooth_start[ida] = this->start(ida) - interp->ndetail_smooth_extend_front();
                    smooth_end[ida]   = this->end(ida) + interp->ndetail_smooth_extend_back();
                }
            }
            // m_log("ibidule = %d -> maks = 1.0 from %d %d %d to %d %d %d", ibidule,
            //       smooth_start[0], smooth_start[1], smooth_start[2],
            //       smooth_end[0], smooth_end[1], smooth_end[2]);
            // apply it
            for_loop(&set_mask_to_one, smooth_start, smooth_end);
        }
    };

    // for each ghost block, set the mask to 1.0 if needed
    iblock_t block_count = 0;
    for (auto gblock : local_parent_) {
        // if ((status_ngh_[block_count] != M_ADAPT_SAME) && (this->status_level() != M_ADAPT_SAME)) {
        //     m_log("block thinks that his neigbor %d (count=%d) has been modified while the block has been changed", gblock->cum_block_id(), block_count);
        //     m_assert(false, "oouuups");
        // }
        mask_smooth(block_count, gblock->ibidule());
        // update the counter
        ++block_count;
    }
    m_assert(block_count == local_parent_.size(), "the two numbers must match: %d vs %ld", block_count, local_parent_.size());
    for (auto gblock : ghost_parent_) {
        // if ((status_ngh_[block_count] != M_ADAPT_SAME) && (this->status_level() != M_ADAPT_SAME)) {
        //     m_log("block thinks that his neigbor %d (count=%d) has been modified while the block has been changed", gblock->cum_block_id(), block_count);
        //     m_assert(false, "oouuups");
        // }
        mask_smooth(block_count, gblock->ibidule());
        // update the counter
        ++block_count;
    }
    m_assert(block_count == (local_parent_.size() + ghost_parent_.size()), "the two numbers must match: %d vs %ld", block_count, (local_parent_.size() + ghost_parent_.size()));

    //................................................
    // smooth depending on the mask
    // SubBlock block_src(this->gs(), this->stride(), -interp->nghost_front(), M_N + interp->nghost_back());
    // SubBlock block_det(this->gs(), this->stride(), -interp->ndetail_smooth_extend_front(), M_N + interp->ndetail_smooth_extend_back());
    SubBlock block_src( -interp->nghost_front(), M_N + interp->nghost_back());
    SubBlock block_det( -interp->ndetail_smooth_extend_front(), M_N + interp->ndetail_smooth_extend_back());

    // do it for every field
    for (auto fid = field_start; fid != field_end; ++fid) {
        auto current_field = fid->second;
        if (!current_field->is_temp()) {
            for (lda_t ida = 0; ida < current_field->lda(); ida++) {
                interp->SmoothOnMask(&block_src, this, this->data(current_field, ida), &block_det, mask);
            }
        }
    }

    // free the status array
    m_free(status_ngh_);
    //-------------------------------------------------------------------------
}

/**
 * @brief clean the points affected by the resolution jump 
 * 
 * @param interp 
 * @param field_start 
 * @param field_end 
 * @param profiler 
 */
void GridBlock::ClearResolutionJump(const Wavelet*  interp, std::map<std::string, Field*  >::const_iterator field_start, std::map<std::string, Field*  >::const_iterator field_end, Prof*  profiler) {
    // //-------------------------------------------------------------------------
    // // reset the temp memory to 0.0
    // memset(coarse_ptr_(), 0, CartBlockMemNum(1) * sizeof(real_t));
    // data_ptr mask      = coarse_ptr_(0, this);
    // real_t*  mask_data = mask.Write();

    // //................................................
    // // lambda to obtain the cleaning pattern
    // auto mask_smooth = [=](const iface_t ibidule) -> void {
    //     // create the lambda to put 1.0
    //     auto set_mask_to_one = [=, &mask_data](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
    //         mask_data[m_idx(i0, i1, i2, 0, this->stride())] = 1.0;
    //     };

    //     // get the sign of the ibidule
    //     real_t sign[3];
    //     GhostGetSign(ibidule, sign);

    //     // create the start and send indexes
    //     bidx_t smooth_start[3], smooth_end[3];
    //     for (lda_t ida = 0; ida < 3; ++ida) {
    //         if (sign[ida] > 0.5) {
    //             // my ngh assumed 0 details in my block
    //             smooth_start[ida] = this->end(ida) - interp->ndetail_citerion_extend_front();
    //             // the number of my ngh details influencing my values
    //             smooth_end[ida] = this->end(ida);
    //         } else if (sign[ida] < (-0.5)) {
    //             // my ngh assumed 0 details in my block
    //             smooth_start[ida] = this->start(ida);
    //             // the number of my ngh details influencing my values
    //             smooth_end[ida] = this->start(ida) + interp->ndetail_citerion_extend_back();
    //         } else {
    //             // even in the directions orthogonal to ibidule, the details must be killed!
    //             // as my neighbor, which might be fine will kill them as well
    //             smooth_start[ida] = this->start(ida);
    //             smooth_end[ida]   = this->end(ida);
    //         }
    //     }
    //     // apply it
    //     for_loop(&set_mask_to_one, smooth_start, smooth_end);
    // };

    // // for each ghost block, set the mask to 1.0 if needed
    // iblock_t block_count = 0;
    // for (auto gblock : local_parent_) {
    //     mask_smooth(gblock->ibidule());
    // }
    // m_assert(block_count == local_parent_.size(), "the two numbers must match: %d vs %ld", block_count, local_parent_.size());
    // for (auto gblock : ghost_parent_) {
    //     mask_smooth(gblock->ibidule());
    // }

    // //................................................
    // // smooth depending on the mask
    // SubBlock block_src(this->gs(), this->stride(), -interp->nghost_front(), M_N + interp->nghost_back());
    // SubBlock block_det(this->gs(), this->stride(), 0, M_N);

    // // do it for every field
    // for (auto fid = field_start; fid != field_end; ++fid) {
    //     auto current_field = fid->second;
    //     for (lda_t ida = 0; ida < current_field->lda(); ida++) {
    //         interp->ClearOnMask(&block_src, this, this->data(current_field, ida), &block_det, mask);
    //     }
    // }

    // // free the status array
    // m_free(status_ngh_);
    // //-------------------------------------------------------------------------
}

/**
 * @brief Compute the detail coefficient and store them in the field details
 * 
 * @param interp the wavelet object
 * @param criterion the criterion field
 * @param details the detail field with the compute detail values
 */
void GridBlock::StoreDetails(const Wavelet*  interp, const Field*  criterion, const Field*  details) {
    m_assert(criterion->ghost_status(), "the field <%s> must have up-to-date ghosts", criterion->name().c_str());
    m_assert(criterion->lda() == details->lda(), "field <%s> and <%s> must have the same size", criterion->name().c_str(), details->name().c_str());
    //-------------------------------------------------------------------------
    for (lda_t ida = 0; ida < criterion->lda(); ida++) {
        SubBlock block_src( -interp->nghost_front(), M_N + interp->nghost_back());
        interp->StoreDetails(&block_src, this->data(criterion, ida), this, this->data(details, ida));
    }
    //-------------------------------------------------------------------------
}

/**
 * @brief resolve the dependency list created while adapting the mesh (see cback_UpdateDependency() ) by interpolating the needed blocks
 * 
 * @param interp the Wavelet object to use for the interpolation/refinement
 * @param field_start the first field to take into account
 * @param field_end the last field to take into account
 */
void GridBlock::SolveDependency(const Wavelet*  interp, std::map<std::string, Field*  >::const_iterator field_start, std::map<std::string, Field*  >::const_iterator field_end, Prof*  profiler) {
    m_assert(n_dependency_active_ == 0 || n_dependency_active_ == 1 || n_dependency_active_ == P8EST_CHILDREN, "wrong value for n_dependency_active_");
    m_profStart(profiler, "solve dependency");
    //-------------------------------------------------------------------------
    if (n_dependency_active_ == 1) {  // this is REFINEMENT
        // if I get only one dependency, I am a child and I need refinement from my parent
        GridBlock* root = this->PopDependency(0);
        m_assert(n_dependency_active_ == 0, "I should be empty now");

        // check the status
        m_assert(this->status_level() == M_ADAPT_NEW_FINE, "my status must be M_ADAPT_NEW_FINE instead of %d", this->status_level());
        m_assert(root->status_level() == M_ADAPT_FINER, "the status of the new root must be M_ADAPT_FINER instead of %d", root->status_level());

        // from the parent, we interpolate to me
        int childid = p4est_GetChildID(xyz_, level_);

        // get the shift given the child id
        const lid_t shift[3]     = {M_NCENTER * ((childid % 2)), M_NCENTER * ((childid % 4) / 2), M_NCENTER * ((childid / 4))};
        const lid_t src_start[3] = {shift[0] - M_GS, shift[1] - M_GS, shift[2] - M_GS};
        const lid_t src_end[3]   = {shift[0] + M_NCENTER + M_GS, shift[1] + M_NCENTER + M_GS, shift[2] + M_NCENTER + M_GS};
        // SubBlock    mem_src(M_GS, M_STRIDE, src_start, src_end);
        SubBlock    mem_src(src_start, src_end);

        // for every field on my parent, interpolate it
        m_assert(mem_map_.size() == 0, "the block should be empty here");
        for (auto fid = field_start; fid != field_end; ++fid) {
            auto current_field = fid->second;
            // allocate the field on me (has not been done yet)
            this->AddField(current_field);
            // refine for every dimension, if the field is not temp
            if (!current_field->is_temp()) {
                for (sid_t ida = 0; ida < current_field->lda(); ida++) {
                    // get the pointers
                    interp->Interpolate(-1, shift, &mem_src, root->data(current_field, ida), this, this->data(current_field, ida));
                }
            }
        }
        // set the status of my newblock
        // remove my ref from the parent
        GridBlock* this_should_be_me = root->PopDependency(childid);
        m_assert(this_should_be_me == this, "this should be me");

        // destroy my parent if I was the last one
        if (root->n_dependency_active() == 0) {
            delete (root);
        }
    } else if (n_dependency_active_ == P8EST_CHILDREN) {  // this is COARSENING
        m_assert(this->status_level() == M_ADAPT_NEW_COARSE, "my status must be M_ADAPT_COARSER instead of %d", this->status_level());
        // I have 8 deps, I am a root, waiting data from coarsening of my children
        //allocate the new fields
        m_assert(mem_map_.size() == 0, "the block should be empty here");
        for (auto fid = field_start; fid != field_end; ++fid) {
            auto current_field = fid->second;
            // allocate the field on me (has not been done yet)
            this->AddField(current_field);
        }

        // I am a parent and I need to fillout my children
        for (sid_t childid = 0; childid < P8EST_CHILDREN; ++childid) {
            GridBlock* child_block = this->PopDependency(childid);
            m_assert(child_block->level() - this->level() == 1, "the child block is not a child");
            m_assert(childid == p4est_GetChildID(child_block->xyz(), child_block->level()), "the two ids must match");
            m_assert(child_block->status_level() == M_ADAPT_COARSER, "the status of the new root must be M_ADAPT_NEW_COARSE instead of %d", child_block->status_level());

            // get the shift for me
            const lid_t shift[3]     = {-M_N * ((childid % 2)), -M_N * ((childid % 4) / 2), -M_N * ((childid / 4))};
            const lid_t trg_start[3] = {M_NCENTER * ((childid % 2)), M_NCENTER * ((childid % 4) / 2), M_NCENTER * ((childid / 4))};
            const lid_t trg_end[3]   = {trg_start[0] + M_NCENTER, trg_start[1] + M_NCENTER, trg_start[2] + M_NCENTER};
            SubBlock    mem_trg(trg_start, trg_end);
            // and an extended source block for my child
            const lid_t src_start[3] = {-M_GS, -M_GS, -M_GS};
            const lid_t src_end[3]   = {M_N + M_GS, M_N + M_GS, M_N + M_GS};
            SubBlock    mem_src(src_start, src_end);

            // for every field, we interpolate it
            for (auto fid = field_start; fid != field_end; ++fid) {
                auto current_field = fid->second;
                if (!current_field->is_temp()) {
                    // interpolate for every dimension
                    for (sid_t ida = 0; ida < current_field->lda(); ida++) {
                        // get the pointers
                        interp->Interpolate(1, shift, &mem_src, child_block->data(current_field, ida), &mem_trg, this->data(current_field, ida));
                    }
                }
            }
            // remove the ref to the child and my ref in the child
            child_block->PopDependency(0);
            m_assert(child_block->n_dependency_active() == 0, "the child must be empty now");
            // delete the block as we are the only one to access it
            delete (child_block);
        }
    } else {
        // we don't refine or coarsen so we compute the Inverse Wavelet Transform
    }
    m_profStop(profiler, "solve dependency");
    //-------------------------------------------------------------------------
}

/**
 * @brief remove the dependency for the given child
 * 
 * @param child_id the id of the dependency to remove, 0 if the parent is removed
 * @return GridBlock* 
 */
GridBlock* GridBlock::PopDependency(const sid_t child_id) {
    m_assert(0 <= child_id && child_id < P8EST_CHILDREN, "child id is out of bound");
    m_assert(dependency_[child_id] != nullptr, "there is nobody here: dep[%d] = %p", child_id, dependency_[child_id]);
    //-------------------------------------------------------------------------
    --n_dependency_active_;
    GridBlock* block      = dependency_[child_id];
    dependency_[child_id] = nullptr;
    return block;
    //-------------------------------------------------------------------------
}

/**
 * @brief add a dependency for a given child
 * 
 * @param child_id the id of the dependency to add, 0 if the dependent is the parent
 * @param dependent_block the pointer to the dependent block
 */
void GridBlock::PushDependency(const sid_t child_id, GridBlock* dependent_block) {
    m_assert(0 <= child_id && child_id < P8EST_CHILDREN, "child id is out of bound");
    m_assert(dependency_[child_id] == nullptr, "there is already someone here");
    //-------------------------------------------------------------------------
    ++n_dependency_active_;
    dependency_[child_id] = dependent_block;
    //-------------------------------------------------------------------------
}

/**
 * @brief build the list of ghosts
 * 
 * @param qid the current quadrant ID
 * @param grid the grid to use to recover the ghosts etc
 * @param interp the wavelet to use (for coarse ghost sizes etc)
 * @param local2disp_window the displacement information for RMA
 */
void GridBlock::GhostInitLists(const qid_t*  qid, const ForestGrid*  grid, const Wavelet*  interp, MPI_Win local2disp_window) {
    //-------------------------------------------------------------------------
    // allocate the ghost pointer, which is reused for the wavelets smoothing
    // size_t alloc_size = m_max(interp->CoarseSize(), m_blockmemsize(1));
    // AllocateCoarsePtr(alloc_size);
    // m_log("I allocate %ld doubles",alloc_size);
    m_assert(interp->CoarseSize() <= CartBlockMemNum(1), "the coarse size must be smaller than a blockmemsize to fit in the coarse memory");
    // m_log("Coarse = %ld vs block size = %ld", interp->CoarseSize(), CartBlockMemNum(1));

    //................................................
    p8est_t*              forest  = grid->p4est_forest();
    p8est_mesh_t*         mesh    = grid->p4est_mesh();
    p8est_ghost_t*        ghost   = grid->p4est_ghost();
    p8est_connectivity_t* connect = forest->connectivity;

    std::list<qdrt_t*>  ngh_list;
    std::list<iblock_t> bid_list;
    std::list<rank_t>   rank_list;

    //................................................
    // get the number of ghost and the min/max of a block
    lid_t  block_min[3], block_max[3];
    real_t block_len[3];
    real_t coarse_hgrid[3];
    for (lda_t id = 0; id < 3; id++) {
        m_assert(level() >= 0, "the level=%d must be >=0", level());
        // set the number of ghost to compute
        block_min[id]    = -interp->nghost_front();
        block_max[id]    = M_N + interp->nghost_back();
        block_len[id]    = p4est_QuadLen(level());
        coarse_hgrid[id] = CoarseHGrid(p4est_QuadLen(level()));
    }

    rank_t my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    // we do the loop in the opposite way, starting with the corners, edges and finally the
    for (iface_t ibidule = (M_NNEIGHBORS - 1); ibidule >= 0; ibidule--) {
        // set the current status to none, whatever happends next
        // ngh_status[ibidule] = NS_NONE;

        //................................................
        p4est_GetNeighbor(forest, connect, ghost, mesh, qid->tid, qid->qid, ibidule, &ngh_list, &bid_list, &rank_list);
        const iblock_t nghosts = ngh_list.size();

        //................................................
        // no ghosts? then is a physical BC
        if (nghosts == 0) {
            // we only apply the physics to entire faces
            if (ibidule < 6) {
                // register the status
                // ngh_status[ibidule] = NS_PHYS;
                PhysBlock* pb = new PhysBlock(ibidule, this, interp->nghost_front(), interp->nghost_back());
                //#pragma omp critical
                phys_.push_back(pb);
                m_verb("I found a physical boundary ghost!\n");
            }
            // else, the edges and corners will be filled through the face
        }

        //................................................
        // this is a real block or a ghost
        for (iblock_t nid = 0; nid < nghosts; ++nid) {
            const qdrt_t*  nghq       = ngh_list.back();
            const iblock_t ngh_cum_id = bid_list.back();
            const rank_t   ngh_rank   = rank_list.back();
            const bool     isghost    = (ngh_rank != my_rank);
            m_verb("reading th list: adress: %p  and rank %d -> is ghost? %d", nghq, ngh_rank, isghost);

            // get the sign, i.e. the normal to the face, the edge of the corner we consider
            real_t sign[3];
            GhostGetSign(ibidule, sign);

            //................................................
            // get the position of the neighbor, as seen by me!!! may be different than the actual one if there is a periodic bc
            real_t ngh_pos[3];
            if (!isghost) {
                // cannot use the p8est function because the which_tree is not assigned, so we retrieve the position through the block
                GridBlock* ngh_block = *(reinterpret_cast<GridBlock**>(nghq->p.user_data));
                ngh_pos[0]           = ngh_block->xyz(0);
                ngh_pos[1]           = ngh_block->xyz(1);
                ngh_pos[2]           = ngh_block->xyz(2);
            } else {
                p8est_qcoord_to_vertex(grid->p4est_connect(), nghq->p.piggy3.which_tree, nghq->x, nghq->y, nghq->z, ngh_pos);
            }
            // fix the shift in coordinates needed if the domain is periodic
            for (lda_t id = 0; id < 3; ++id) {
                // if we are periodic, we overwrite the position in the direction of the normal !!ONLY!!
                // since it is my neighbor in this normal direction, I am 100% sure that it's origin corresponds to the end of my block
                const real_t to_replace = sign[id] * sign[id] * grid->domain_periodic(id);  // is (+-1)^2 = +1 if we need to replace it, 0.0 otherwize
                // get the expected position
                m_assert(level() >= 0, "the level=%d must be >=0", level());
                m_assert(nghq->level >= 0, "the level=%d must be >=0", nghq->level);
                const real_t expected_pos = xyz(id) + (sign[id] > 0.5) * p4est_QuadLen(level()) - (sign[id] < -0.5) * p4est_QuadLen(nghq->level);
                // we override the position if a replacement is needed only
                ngh_pos[id] = to_replace * expected_pos + (1.0 - to_replace) * ngh_pos[id];
            }
            // get the hgrid
            m_assert(nghq->level >= 0, "the level=%d must be >=0", nghq->level);
            const real_t ngh_len[3] = {p4est_QuadLen(nghq->level), p4est_QuadLen(nghq->level), p4est_QuadLen(nghq->level)};
            // const real_t ngh_hgrid[3] = {p4est_QuadLen(nghq->level) / (M_N-1), p4est_QuadLen(nghq->level) / M_N, p4est_QuadLen(nghq->level) / M_N};
            const real_t ngh_hgrid[3] = {CartBlockHGrid(ngh_len[0]), CartBlockHGrid(ngh_len[1]), CartBlockHGrid(ngh_len[2])};

            //................................................
            // create the new block and push back
            if (!isghost) {
                // m_log("the block is not a ghost");
                // associate the corresponding neighboring block
                GridBlock* ngh_block = p4est_GetGridBlock(nghq);
                // #ifndef NDEBUG
                // {  // check the indexing... you never know with this shitty functions
                //     m_log("check the indexing...");
                //     m_log("try to get tree number %d", nghq->p.piggy3.which_tree);
                //     p8est_tree_t*  tree    = p8est_tree_array_index(forest->trees, nghq->p.piggy1.which_tree);
                //     p4est_locidx_t quad_id = ngh_cum_id - tree->quadrants_offset;
                //     m_assert(quad_id >= 0, "the quad id must be >0: %d = %d - %d", ngh_cum_id, tree->quadrants_offset);
                //     p8est_quadrant* quad = p8est_quadrant_array_index(&tree->quadrants, quad_id);
                //     m_assert(p4est_GetGridBlock(quad) == ngh_block, "these two addresses must be the same! %p vs %p", p4est_GetGridBlock(quad), ngh_block);
                //     m_log("end of check, compute the intersections");
                // }
                // #endif
                // m_log("block @ %f %f %f: neighbor num %d @ level = %d", xyz_[0], xyz_[1], xyz_[2], ibidule, ngh_block->level());

                // register the gb in a list
                if (nghq->level == level()) {
                    // m_log("creating a same level");
                    // sibling: source = neighbor GridBlock, target = me
                    GBLocal* gb = new GBLocal(ibidule, ngh_cum_id);
                    gb->Intersect(/* source */ ngh_block->level(), ngh_pos, ngh_hgrid, ngh_len,
                                  /* traget */ level_, xyz_, hgrid_, block_min, block_max);
                    gb->data_src(ngh_block);

                    //#pragma omp critical
                    local_sibling_.push_back(gb);

                } else if (nghq->level < level()) {
                    // m_log("creating a coarser");
                    // parent: source = neighbor, target = me
                    GBLocal* gb = new GBLocal(ibidule, ngh_cum_id);
                    gb->Intersect(/* source */ ngh_block->level(), ngh_pos, ngh_hgrid, ngh_len,
                                  /* target */ level_, xyz_, hgrid_, block_min, block_max);
                    gb->data_src(ngh_block);
                    //#pragma omp critical
                    local_parent_.push_back(gb);

                    // m_log("creating a reverse");
                    // the children: the source = the coarse myself, target = my neighbor
                    GBLocal* invert_gb = new GBLocal(-1, ngh_cum_id);
                    invert_gb->Intersect(/* source */ level_ - 1, xyz_, coarse_hgrid, block_len,
                                         /* target */ ngh_block->level(), ngh_pos, ngh_hgrid, block_min, block_max);
                    invert_gb->data_src(ngh_block);
                    //#pragma omp critical
                    local_parent_reverse_.push_back(invert_gb);
                } else if (nghq->level > level()) {
                    // m_log("creating a finer");
                    m_assert((nghq->level - level_) == 1, "The delta level is not correct: %d - %d", nghq->level, level_);
                    // register the coarse
                    GBLocal* gb = new GBLocal(ibidule, ngh_cum_id);
                    //#pragma omp critical
                    local_children_.push_back(gb);
                } else {
                    m_assert(false,"this shouldn't happen");
                }
            }
            //................................................
            else {
                // get the local number in the remote rank and the remote rank
                // iblock_t ngh_local_id = nghq->p.piggy3.local_num;
                m_assert(ngh_cum_id == nghq->p.piggy3.local_num, "the two numbering must match: %d vs %d", ngh_cum_id, nghq->p.piggy3.local_num);
                // rank_t ngh_rank     = p4est_GetOwnerFromGhost(forest, nghq);
                m_assert(ngh_rank >= 0, "p4est unable to recover the rank... baaaad news: %d", ngh_rank);
                m_assert(ngh_rank < forest->mpisize, "the rank must be smaller than the comm size: %d vs %d ", ngh_rank, forest->mpisize);
                m_assert((forest->global_first_quadrant[ngh_rank + 1] - forest->global_first_quadrant[ngh_rank]) > 0, "the neighbor must have quadrants");
                m_assert(0 <= ngh_cum_id, "the piggy3.local_num must fall in the quad's limits: %d < %d", 0, ngh_cum_id);
                m_assert(ngh_cum_id < (forest->global_first_quadrant[ngh_rank + 1] - forest->global_first_quadrant[ngh_rank]), "the piggy3.local_num must fall in the quad's limits: %d < %ld", ngh_cum_id, (forest->global_first_quadrant[ngh_rank + 1] - forest->global_first_quadrant[ngh_rank]));

                // register the ghost block in a list
                //................................................
                if (nghq->level == level_) {
                    // sibling: source = neighbor GridBlock, target = me
                    GBMirror* gb = new GBMirror(ibidule, ngh_cum_id, ngh_rank);
                    gb->Intersect(/* source */ nghq->level, ngh_pos, ngh_hgrid, ngh_len,
                                  /* target */ level(), xyz(), hgrid(), block_min, block_max);
                    // ask the displacement (will be available later, when completing the call)
                    m_assert(ngh_cum_id >= 0, "the ngh_local_id = %d must be non-negative", ngh_cum_id);
                    MPI_Get(gb->data_src_ptr(), 1, MPI_AINT, ngh_rank, ngh_cum_id, 1, MPI_AINT, local2disp_window);
                    //#pragma omp critical
                    ghost_sibling_.push_back(gb);

                }
                //................................................
                else if (nghq->level < level_) {
                    // parent: source = neighbor, target = me
                    GBMirror* gb = new GBMirror(ibidule, ngh_cum_id, ngh_rank);
                    gb->Intersect(/* source */ nghq->level, ngh_pos, ngh_hgrid, ngh_len,
                                  /* target */ level(), xyz(), hgrid(), block_min, block_max);
                    // ask the displacement (will be available later, when completing the call
                    m_assert(ngh_cum_id >= 0, "the ngh_local_id = %d must be non-negative", ngh_cum_id);
                    MPI_Get(gb->data_src_ptr(), 1, MPI_AINT, ngh_rank, ngh_cum_id, 1, MPI_AINT, local2disp_window);
                    //#pragma omp critical
                    ghost_parent_.push_back(gb);

                    // I compute my own contribution to my neighbor ghost points
                    GBMirror* invert_gb = new GBMirror(ibidule, ngh_cum_id, ngh_rank);
                    invert_gb->Intersect(/* source */ level() - 1, xyz(), coarse_hgrid, block_len,
                                         /* target */ nghq->level, ngh_pos, ngh_hgrid, block_min, block_max);
                    m_assert(ngh_cum_id >= 0, "the ngh_local_id = %d must be non-negative", ngh_cum_id);
                    MPI_Get(invert_gb->data_src_ptr(), 1, MPI_AINT, ngh_rank, ngh_cum_id, 1, MPI_AINT, local2disp_window);
                    //#pragma omp critical
                    ghost_parent_reverse_.push_back(invert_gb);

                }
                //................................................
                else if (nghq->level > level_) {
                    const real_t ngh_hgrid_coarse[3] = {CoarseHGrid(ngh_len[0]), CoarseHGrid(ngh_len[1]), CoarseHGrid(ngh_len[2])};

                    // children: source = coarse version of my neighbor, target = myself
                    GBMirror* gb = new GBMirror(ibidule, ngh_cum_id, ngh_rank);
                    gb->Intersect(/* source */ nghq->level - 1, ngh_pos, ngh_hgrid_coarse, ngh_len,
                                  /* target */ level(), xyz(), hgrid(), block_min, block_max);
                    //#pragma omp critical
                    ghost_children_.push_back(gb);
                }
                //................................................
                else {
                    m_assert(false, "The delta level is not correct: %d - %d", nghq->level, level());
                }
            }
            // pop the last element, will go to the next one
            ngh_list.pop_back();
            rank_list.pop_back();
            bid_list.pop_back();
        }
    }
    // //-------------------------------------------------------------------------
    // for (iface_t ibidule = (M_NNEIGHBORS - 1); ibidule >= 0; ibidule--) {
    //     // set the current status to
    //     ngh_status[ibidule] = NS_NONE;
    // }
}

/**
 * @brief free and clear the ghost list
 * 
 */
void GridBlock::GhostFreeLists() {
    //-------------------------------------------------------------------------
    // need to free the memory as well, otherwise I cannot reallocate it later on
    // coarse_ptr_.Free();

    // clear the ghost lists
    auto remove_block = [](auto block) { delete (block); };
    std::for_each(local_sibling_.begin(), local_sibling_.end(), remove_block);
    std::for_each(local_parent_.begin(), local_parent_.end(), remove_block);
    std::for_each(local_children_.begin(), local_children_.end(), remove_block);
    std::for_each(local_parent_reverse_.begin(), local_parent_reverse_.end(), remove_block);
    std::for_each(ghost_sibling_.begin(), ghost_sibling_.end(), remove_block);
    std::for_each(ghost_parent_.begin(), ghost_parent_.end(), remove_block);
    std::for_each(ghost_children_.begin(), ghost_children_.end(), remove_block);
    std::for_each(ghost_parent_reverse_.begin(), ghost_parent_reverse_.end(), remove_block);
    std::for_each(phys_.begin(), phys_.end(), remove_block);

    // clear the lists
    local_sibling_.clear();
    local_parent_.clear();
    local_children_.clear();
    local_parent_reverse_.clear();
    ghost_sibling_.clear();
    ghost_parent_.clear();
    ghost_children_.clear();
    ghost_parent_reverse_.clear();
    phys_.clear();
    //-------------------------------------------------------------------------
}

/**
 * @brief Do the first part of GhostGet, everything which can be done without MPI
 * 
 * @warning we assume that GhostGet_Post has been called first
 * 
 * - get the siblings values
 * - get the coarser values to my coarse temporary array
 * 
 * @param field the field to interpolate
 * @param ida the current dimension
 * @param interp the wavelet
 * @param mirrors_window the window where to find the mirrors
 */
void GridBlock::GhostGet_Cmpt(const Field*  field, const lda_t ida, const Wavelet*  interp) {
    //-------------------------------------------------------------------------
    // get the siblings
    {
        const SubBlock bsrc_neighbor(0, M_N);
        const data_ptr data_trg = data(field, ida);

        // get the copy of local siblings
        for (const auto gblock : (local_sibling_)) {
            GridBlock*       ngh_block = gblock->data_src();
            const_data_ptr   data_src  = const_data_ptr(ngh_block->data(field, ida));
            const Layout* block_trg = gblock;
            // copy the information
            interp->Copy(gblock->dlvl(), gblock->shift(), &bsrc_neighbor, data_src, block_trg, data_trg);
        }
    }

    //................................................
    const bool do_coarse = (local_parent_.size() + ghost_parent_.size()) > 0;
    if (do_coarse) {
        // we use the full neighbor span
        const SubBlock bsrc_neighbor(0, M_N);

        // copy the siblings to coarse
        for (const auto gblock : local_sibling_) {
            GridBlock* ngh_block = gblock->data_src();
            SubBlock   block_trg;
            interp->CoarseFromFine(gblock, &block_trg);
            const_data_ptr data_src = const_data_ptr(ngh_block->data(field, ida));
            data_ptr data_trg = coarse_ptr_(0, interp->CoarseNGhostFront(),interp->CoarseStride());
            // interpolate, the level is 1 coarser and the shift is unchanged
            m_assert((gblock->dlvl() + 1) == 1, "the difference of level MUST be 1");
            interp->Copy(gblock->dlvl() + 1, gblock->shift(), &bsrc_neighbor, data_src, &block_trg, data_trg);
        }
        // copy the parents to coarse
        for (const auto gblock : local_parent_) {
            GridBlock* ngh_block = gblock->data_src();
            SubBlock   block_trg;
            interp->CoarseFromFine(gblock, &block_trg);
            const_data_ptr data_src = const_data_ptr(ngh_block->data(field, ida));
            data_ptr data_trg = coarse_ptr_(0, interp->CoarseNGhostFront(),interp->CoarseStride());
            // interpolate, the level is 1 coarser and the shift is unchanged
            m_assert((gblock->dlvl() + 1) == 0, "the difference of level MUST be 1");
            interp->Copy(gblock->dlvl() + 1, gblock->shift(), &bsrc_neighbor, data_src, &block_trg, data_trg);
        }
    }
    //-------------------------------------------------------------------------
}

/**
 * @brief Do the first part of GhostGet, everything that requires MPI
 * 
 * @param field 
 * @param ida 
 * @param interp 
 * @param mirrors_window 
 */
void GridBlock::GhostGet_Post(const Field*  field, const lda_t ida, const Wavelet*  interp, MPI_Win mirrors_window) {
    //-------------------------------------------------------------------------
    // get the siblings
    {
        const_data_ptr bsrc_neighbor(nullptr, M_GS, M_STRIDE);
        // const SubBlock bsrc_neighbor(0, M_N);
        const data_ptr data_trg = data(field, ida);

        // start the RMA for the sibling neighbors
        for (const auto gblock : ghost_sibling_) {
            const MPI_Aint   disp_src  = gblock->data_src();
            const rank_t     disp_rank = gblock->rank();
            const Layout* block_trg = gblock;
            // copy the information
            interp->GetRma(gblock->dlvl(), gblock->shift(), bsrc_neighbor, disp_src, block_trg, data_trg, disp_rank, mirrors_window);
        }
    }
    //................................................
    const bool do_coarse = (local_parent_.size() + ghost_parent_.size()) > 0;
    if (do_coarse) {
        // reset the coarse memory
        memset(coarse_ptr_(), 0, interp->CoarseSize() * sizeof(real_t));
        // we use the full neighbor span
        const_data_ptr bsrc_neighbor(nullptr, M_GS, M_STRIDE);

        //................................................
        // RMA the sibligns ghosts to the tmp
        for (const auto gblock : ghost_sibling_) {
            MPI_Aint disp_src  = gblock->data_src();
            rank_t   disp_rank = gblock->rank();
            // get the associated coarse block
            SubBlock block_trg;
            interp->CoarseFromFine(gblock, &block_trg);
            data_ptr data_trg = coarse_ptr_(0, interp->CoarseNGhostFront(), interp->CoarseStride());  // + m_zeroidx(0, &block_trg);
            // interpolate, the level is 1 coarser and the shift is unchanged
            m_assert((gblock->dlvl() + 1) == 1, "the difference of level MUST be 1 or 0");
            interp->GetRma((gblock->dlvl() + 1), gblock->shift(), bsrc_neighbor, disp_src, &block_trg, data_trg, disp_rank, mirrors_window);
        }
        // RMA the parent ghosts to the tmp
        for (const auto gblock : ghost_parent_) {
            MPI_Aint disp_src  = gblock->data_src();
            rank_t   disp_rank = gblock->rank();
            // get the associated coarse block
            SubBlock block_trg;
            interp->CoarseFromFine(gblock, &block_trg);
            data_ptr data_trg = coarse_ptr_(0, interp->CoarseNGhostFront(), interp->CoarseStride());
            // interpolate, the level is 1 coarser and the shift is unchanged
            m_assert((gblock->dlvl() + 1) == 0, "the difference of level MUST be 1 or 0");
            interp->GetRma((gblock->dlvl() + 1), gblock->shift(), bsrc_neighbor, disp_src, &block_trg, data_trg, disp_rank, mirrors_window);
        }
    }
    //-------------------------------------------------------------------------
}

/**
 * @brief Do the second part of the GhostGet
 * 
 * @warning First call @ref GhostGet_Post() and @ref GhostGet_Cmpt()
 * 
 * - wait for the communications to be over
 * - add the physical BC to my coarse temp
 * - interpolate to my ghost points
 * 
 * @param field 
 * @param ida 
 * @param interp 
 */
void GridBlock::GhostGet_Wait(const Field*  field, const lda_t ida, const Wavelet*  interp) {
    //-------------------------------------------------------------------------
    const bool do_coarse = (local_parent_.size() + ghost_parent_.size()) > 0;
    if (do_coarse) {
        //................................................
        {
            // m_log("reset the coarse block to %d and stride %d ", interp->CoarseNGhostFront(), interp->CoarseStride());
            const SubBlock coarse_block(0, M_NHALF);
            // copy myself to the coarse, one point out of 2
            const lid_t    shift[3] = {0, 0, 0};
            const_data_ptr data_src = const_data_ptr(data(field, ida));
            const data_ptr data_trg = coarse_ptr_(0, interp->CoarseNGhostFront(), interp->CoarseStride());
            // interpolate
            interp->Copy(1, shift, this, data_src, &coarse_block, data_trg);
        }

        //................................................
        // do here some physics, to completely fill the coarse block before the interpolation
        for (auto gblock : phys_) {
            // get the direction and the corresponding bctype
            const bctype_t bctype = field->bctype(ida, gblock->iface());
            // in the face direction, the start and the end are already correct, only the fstart changes
            SubBlock coarse_block;
            interp->CoarseFromFine(gblock, &coarse_block);
            lid_t fstart[3];
            interp->CoarseFromFine(face_start[gblock->iface()], fstart);
            data_ptr data_trg = coarse_ptr_(0, interp->CoarseNGhostFront(), interp->CoarseStride());
            // apply the BC
            // m_log("apply bc on coarse for block @ %f %f %f", xyz(0), xyz(1), xyz(2));
            if (bctype == M_BC_NEU) {
                NeumanBoundary<M_WAVELET_N - 1> bc;
                bc(gblock->iface(), fstart, hgrid_, 0.0, &coarse_block, data_trg);
            } else if (bctype == M_BC_DIR) {
                DirichletBoundary<M_WAVELET_N - 1> bc;
                bc(gblock->iface(), fstart, hgrid_, 0.0, &coarse_block, data_trg);
            } else if (bctype == M_BC_EXTRAP) {
                ExtrapBoundary<M_WAVELET_N> bc;
                bc(gblock->iface(), fstart, hgrid_, 0.0, &coarse_block, data_trg);
            } else if (bctype == M_BC_ZERO) {
                ZeroBoundary bc;
                bc(gblock->iface(), fstart, hgrid_, 0.0, &coarse_block, data_trg);
            } else {
                m_assert(false, "this type of BC is not implemented yet or not valid %d for field <%s>", bctype, field->name().c_str());
            }
        }
        //................................................
        // refine from the coarse to the parents
        {

            // take the full coarse block and set the info in my GP
            const SubBlock block_src(-interp->CoarseNGhostFront(), M_NHALF + interp->CoarseNGhostBack());
            const_data_ptr data_src = const_data_ptr(coarse_ptr_(0, interp->CoarseNGhostFront(), interp->CoarseStride()));
            lid_t          shift[3] = {0, 0, 0};

            for (const auto gblock : local_parent_) {
                // extension from the ghost
                // gblock->ExtendGhost(interp->ndetail_citerion_extend_front(), interp->ndetail_citerion_extend_back(), &block_trg);
                // interp->Interpolate(-1, shift, &block_src, data_src, &block_trg, data(field, ida));
                // no extension
                interp->Interpolate(-1, shift, &block_src, data_src, gblock, data(field, ida));
            }
            for (const auto gblock : ghost_parent_) {
                interp->Interpolate(-1, shift, &block_src, data_src, gblock, data(field, ida));
            }
        }
    }
    //-------------------------------------------------------------------------
}

/**
 * @brief Do the first part of GhostPut: start the RMA put communications
 * 
 * - apply the physical BC to the best of my knowledge
 * - coarsen my data to obtain the scaling coefficient (interpolate, not copy!)
 * - put the obtained values to my neighbor's ghost points
 * 
 * @param field 
 * @param ida 
 * @param interp 
 * @param mirrors_window 
 */
void GridBlock::GhostPut_Post(const Field*  field, const lda_t ida, const Wavelet*  interp, MPI_Win mirrors_window) {
    //-------------------------------------------------------------------------
    const bool do_coarse = (local_parent_.size() + ghost_parent_.size()) > 0;
    if (do_coarse) {
        //................................................
        // apply the physics to the best of my knowledge
        // some ghosts are missing but this is okay
        data_ptr data_trg = data(field, ida);

        // m_log("apply bc for block @ %f %f %f", xyz(0), xyz(1), xyz(2));
        for (auto gblock : phys_) {
            bctype_t bctype = field->bctype(ida, gblock->iface());
            // get the correct face_start
            if (bctype == M_BC_NEU) {
                NeumanBoundary<M_WAVELET_N - 1> bc;
                bc(gblock->iface(), face_start[gblock->iface()], hgrid_, 0.0, gblock, data_trg);
            } else if (bctype == M_BC_DIR) {
                DirichletBoundary<M_WAVELET_N - 1> bc;
                bc(gblock->iface(), face_start[gblock->iface()], hgrid_, 0.0, gblock, data_trg);
            } else if (bctype == M_BC_EXTRAP) {
                ExtrapBoundary<M_WAVELET_N> bc;
                bc(gblock->iface(), face_start[gblock->iface()], hgrid_, 0.0, gblock, data_trg);
            } else if (bctype == M_BC_ZERO) {
                ZeroBoundary bc;
                bc(gblock->iface(), face_start[gblock->iface()], hgrid_, 0.0, gblock, data_trg);
            } else {
                m_assert(false, "this type of BC is not implemented yet or not valid %d", bctype);
            }
        }

        //................................................
        // enforce the 0-detail policy in the needed region
        {
            memset(coarse_ptr_(), 0, CartBlockMemNum(1) * sizeof(real_t));
            data_ptr mask      = coarse_ptr_(0);
            real_t*  mask_data = mask.Write();

            // lambda to obtain the detail checking pattern
            auto overwrite = [=](const iface_t ibidule) -> void {
                // get the sign of the ibidule
                real_t sign[3];
                GhostGetSign(ibidule, sign);

                // create the start and end indexes
                bidx_t smooth_start[3], smooth_end[3];
                for (lda_t lda = 0; lda < 3; ++lda) {
                    if (sign[lda] > 0.5) {
                        // my ngh assumed 0 details in my block
                        smooth_start[lda] = this->end(lda) - interp->ndetail_citerion_extend_front();
                        // the number of my ngh details influencing my values
                        smooth_end[lda] = this->end(lda);
                    } else if (sign[lda] < (-0.5)) {
                        // my ngh assumed 0 details in my block
                        smooth_start[lda] = this->start(lda);
                        // the number of my ngh details influencing my values
                        smooth_end[lda] = this->start(lda) + interp->ndetail_citerion_extend_back();
                    } else {
                        // even in the directions orthogonal to ibidule, the details must be killed!
                        // as my neighbor, which might be fine will kill them as well
                        smooth_start[lda] = this->start(lda);
                        smooth_end[lda]   = this->end(lda);
                    }
                }
                // apply it
                // for_loop(&set_mask_to_one, smooth_start, smooth_end);
                // compute the detail, store them in the mask
                SubBlock block_src(-interp->nghost_front(), M_N + interp->nghost_back());
                SubBlock block_trg(smooth_start, smooth_end);
                interp->OverwriteDetails(&block_src, &block_trg, this->data(field, ida));
            };
            // let's god
            for (auto gblock : local_parent_) {
                overwrite(gblock->ibidule());
            }
            for (auto gblock : ghost_parent_) {
                overwrite(gblock->ibidule());
            }
        }
        //................................................
        // coarsen to my finer neighbors
        {
            //................................................
            // reset the tmp to use for the put operations
            // coarse_ptr.Reset(interp->CoarseSize(), interp->CoarseNGhostFront(), interp->CoarseStride());
            memset(coarse_ptr_(), 0, interp->CoarseSize() * sizeof(real_t));
            
            // // I am now complete (except children GP), get my coarse representation
            const SubBlock coarse_block(0, M_NHALF);
            data_ptr       data_coarse = coarse_ptr_(0, interp->CoarseNGhostFront(), interp->CoarseStride());

            // the source block is the ghost extended block
            const lid_t    shift[3] = {0, 0, 0};
            const SubBlock me_extended(-interp->nghost_front(), M_N + interp->nghost_back());
            // interpolate, the level is 1 coarser and the shift is unchanged
            interp->Interpolate(1, shift, &me_extended, data(field, ida), &coarse_block, data_coarse);

            //................................................
            // start the put commands to set the value to my neighbors
            // loop on the ghost list
            for (const auto gblock : ghost_parent_reverse_) {
                MPI_Aint disp_trg = gblock->data_src();
                rank_t   trg_rank = gblock->rank();
                const_data_ptr trg_data(nullptr, this->data(field, ida).gs(),this->data(field, ida).stride() );
                // interpolate, the parent's mirror have been created to act on the tmp
                m_assert(gblock->dlvl() == 0, "we must have a level 0, here %d", gblock->dlvl());
                interp->PutRma(gblock->dlvl(), gblock->shift(), &coarse_block, data_coarse, gblock, trg_data, disp_trg, trg_rank, mirrors_window);
            }
            for (const auto gblock : local_parent_reverse_) {
                GridBlock* ngh_block = gblock->data_src();
                data_ptr   data_trg  = ngh_block->data(field, ida);
                // interpolate, the level is 1 coarser and the shift is unchanged
                m_assert(gblock->dlvl() == 0, "we must have a level 0, here %d", gblock->dlvl());
                interp->Copy(gblock->dlvl(), gblock->shift(), &coarse_block, data_coarse, gblock, data_trg);
            }
        }
    }
    //-------------------------------------------------------------------------
}

/**
 * @brief Do the second part of GhostPut
 * 
 * - the comm are over, I now have everything I need to compute the physical BC
 * 
 * @param field 
 * @param ida 
 * @param interp 
 */
void GridBlock::GhostPut_Wait(const Field*  field, const lda_t ida, const Wavelet*  interp) {
    //-------------------------------------------------------------------------
    data_ptr data_trg = data(field, ida);
    // m_log("apply bc for block @ %f %f %f", xyz(0), xyz(1), xyz(2));
    for (auto gblock : phys_) {
        bctype_t bctype = field->bctype(ida, gblock->iface());
        // get the correct face_start
        if (bctype == M_BC_NEU) {
            NeumanBoundary<M_WAVELET_N - 1> bc;
            bc(gblock->iface(), face_start[gblock->iface()], hgrid_, 0.0, gblock, data_trg);
        } else if (bctype == M_BC_DIR) {
            DirichletBoundary<M_WAVELET_N - 1> bc;
            bc(gblock->iface(), face_start[gblock->iface()], hgrid_, 0.0, gblock, data_trg);
        } else if (bctype == M_BC_EXTRAP) {
            ExtrapBoundary<M_WAVELET_N> bc;
            bc(gblock->iface(), face_start[gblock->iface()], hgrid_, 0.0, gblock, data_trg);
        } else if (bctype == M_BC_ZERO) {
            ZeroBoundary bc;
            bc(gblock->iface(), face_start[gblock->iface()], hgrid_, 0.0, gblock, data_trg);
        } else {
            m_assert(false, "this type of BC is not implemented yet or not valid %d", bctype);
        }
        //-------------------------------------------------------------------------
    }
    //-------------------------------------------------------------------------
}

// /**
//  * @brief Downsample the block and re-evaluate the boundary conditions on the coarse version
//  *
//  * @param field
//  * @param ida
//  * @param interp
//  */
// void GridBlock::Coarse_DownSampleWithBoundary(const Field* field, const lda_t ida, const Wavelet* interp, SubBlock* coarse_block) {
//     //-------------------------------------------------------------------------
//     // reset the tmp value
//     memset(coarse_ptr_, 0, interp->CoarseStride());
//     // get the coarse and extended SubBlocks
//     const lid_t n_ghost_front = (interp->ncriterion_front() + 1) / 2;
//     const lid_t n_ghost_back  = (interp->ncriterion_back() + 1) / 2;  // if we have 3 ghosts, we need 2 coarse and not 1...
//     const lid_t stride        = n_ghost_back + n_ghost_front + M_HN;
//     coarse_block->Reset(n_ghost_front, stride, -n_ghost_front, M_HN + n_ghost_back);
//     m_assert(stride <= interp->CoarseStride(), "ohoh the Coarse block is too small: %d <= %ld", stride, interp->CoarseStride());

//     const SubBlock extended_src(M_GS, M_STRIDE, -M_GS, M_N + M_GS);
//     const lid_t    shift[3] = {0, 0, 0};
//     const data_ptr data_trg = coarse_ptr_ + m_zeroidx(0, coarse_block);
//     const data_ptr data_src = data(field, ida);

//     // interpolate
//     interp->Copy(1, shift, &extended_src, data_src, coarse_block, data_trg);

//     //................................................
//     // over-write the BC
//     for (auto gblock : phys_) {
//         m_assert(false, "we shouldn't be here");
//         // get the direction and the corresponding bctype
//         const bctype_t bctype = field->bctype(ida, gblock->iface());
//         // in the face direction, the start and the end are already correct, only the fstart changes
//         SubBlock coarse_block;
//         interp->CoarseFromFine(gblock, &coarse_block);
//         lid_t fstart[3];
//         interp->CoarseFromFine(face_start[gblock->iface()], fstart);
//         // apply the BC
//         if (bctype == M_BC_NEU) {
//             NeumanBoundary<M_WAVELET_N - 1> bc;
//             bc(gblock->iface(), fstart, hgrid_, 0.0, &coarse_block, data_trg);
//         } else if (bctype == M_BC_DIR) {
//             DirichletBoundary<M_WAVELET_N - 1> bc;
//             bc(gblock->iface(), fstart, hgrid_, 0.0, &coarse_block, data_trg);
//         } else if (bctype == M_BC_EXTRAP) {
//             ExtrapBoundary<M_WAVELET_N> bc;
//             bc(gblock->iface(), fstart, hgrid_, 0.0, &coarse_block, data_trg);
//         } else if (bctype == M_BC_ZERO) {
//             ZeroBoundary bc;
//             bc(gblock->iface(), fstart, hgrid_, 0.0, &coarse_block, data_trg);
//         } else {
//             m_assert(false, "this type of BC is not implemented yet or not valid %d", bctype);
//         }
//     }
//     //-------------------------------------------------------------------------
// }
