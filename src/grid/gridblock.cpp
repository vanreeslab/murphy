#include "gridblock.hpp"

#include <p8est_bits.h>

#include <algorithm>
#include <string>

#include "core/forloop.hpp"
#include "core/macros.hpp"
#include "grid/boundary.hpp"
#include "p8est_iterate.h"
#include "tools/toolsp4est.hpp"
#include "tools/toolsmpi.hpp"
#include "core/memdata.hpp"

using std::string;

static const lid_t face_start[6][3] = {{0, 0, 0}, {M_N, 0, 0}, {0, 0, 0}, {0, M_N, 0}, {0, 0, 0}, {0, 0, M_N}};

/**
 * @brief return the grid spacing of a coarser block
 * 
 * @param len the physical length of the considered block
 */
constexpr real_t CoarseHGrid(const real_t len) {
    return len / (real_t)(M_NHALF);
}

//==============================================================================
/**
 * @brief constructs a new GridBlock
 * 
 * @param length the length of the current block
 * @param xyz the position of the origin, i.e. the left,bottom corner, (x,y,z)
 * @param level the level of the block
 */
GridBlock::GridBlock(const real_t length, const real_t xyz[3], const sid_t level) : CartBlock(length, xyz, level) {
    m_begin;
    //--------------------------------------------------------------------------
    status_lvl_ = M_ADAPT_SAME;

    // init the dependencies
    n_dependency_active_ = 0;
    for (sid_t id = 0; id < P8EST_CHILDREN; ++id) {
        dependency_[id] = nullptr;
    }
    // allocate the coarse ptr
    MemLayout myself = BlockLayout();
    coarse_ptr_.Allocate(myself.n_elem);
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief Destroys the GridBlock
 */
GridBlock::~GridBlock() {
    //-------------------------------------------------------------------------
    coarse_ptr_.Free();
    //-------------------------------------------------------------------------
}

/**
 * @brief compute the wavelet refinement/coarsening criterion and update the block status accordingly
 * 
 * The criterion of the block must always be:
 *  ctol <= criterion <= rtol
 * 
 * - if citerion > rtol, we must refine, whatever the other field's dimension value
 * - if citerion < ctol, we can coarsen if this is satisfy by every dimension of the field
 * 
 * @warning the call to this function is done on a "dimension by dimension" basis on the criterion field. The current dimension is given by ida
 * 
 */
void GridBlock::UpdateStatusFromCriterion(const Wavelet* interp, const real_t rtol, const real_t ctol, const Field* field_citerion, const lda_t ida){
    //--------------------------------------------------------------------------
    // const bidx_t ghost_len_interp[2] = {interp->nghost_front(), interp->nghost_back()};
    m_assert(rtol > ctol, "the refinement tolerance must be > the coarsening tolerance: %e vs %e", rtol, ctol);
    //--------------------------------------------------------------------------
    // prevent coarsening of a block that has been refined in the past
    // -> we can always re-refine a block that has been coarsened
    // -> as a matter of fact, if we compute the details on a block that has been refined, they should be 0, which would not make sense to coarsen again.
    const bool forbid_coarsening = status_refined_;
    const bool forbid_refinement = false;
    // determine if I have already a decision done = no corsening + no refinement possible or other dimension decided to refine
    // if another dimension has decided to coarsen, we can always change our mind
    const bool is_over = (forbid_coarsening && forbid_refinement) || (status_lvl_ == M_ADAPT_FINER);

    // if no decision has been made, go for the computation in the current dimension
    if (!is_over) {
        // go to the computation
        // const SubBlock block_src(this->gs(), this->stride(), -interp->nghost_front(), M_N + interp->nghost_back());
        // const SubBlock block_detail(this->gs(), this->stride(), -interp->ndetail_citerion_extend_front(), M_N + interp->ndetail_citerion_extend_back());
        // const real_t   norm = interp->Criterion(&block_src, this->data(field_citerion, ida), &block_detail);
        const MemSpan      span_src(-interp->nghost_front(), M_N + interp->nghost_back());
        const MemSpan      span_detail(-interp->ndetail_citerion_extend_front(), M_N + interp->ndetail_citerion_extend_back());
        const ConstMemData data_src = this->ConstData(field_citerion, ida);
        const real_t       norm     = interp->Criterion(&span_src, &data_src, &span_detail);
        //(&block_src, this->data(field_citerion, ida), &block_detail);

        // get what we should do = what is safe to do considering this direction
        const bool should_refine  = (norm > rtol) && (!forbid_refinement);
        const bool should_coarsen = (norm < ctol) && (!forbid_coarsening);
        const bool should_stay    = !(should_coarsen || should_refine);
        m_assert((should_coarsen + should_refine + should_stay) == 1, "the sum of the three bools must be 1 = we must make our mind here");

        //......................................................................
        // 1. refinement?
        // if we should refine, we always refine, whatever the other directions have said
        status_lvl_ = (should_refine) ? (M_ADAPT_FINER) : (status_lvl_);
        //......................................................................
        // 2. shouldn't change?
        // if we should stay and the previous directions have said we should coarsen, we cannot coarsen anymore
        status_lvl_ = (should_stay && status_lvl_ == M_ADAPT_COARSER) ? M_ADAPT_SAME : status_lvl_;
        //......................................................................
        // 3. coarsening
        // if we should coarsen and the previous directions said stay the same (status is M_ADAPT_SAME) we coarsen
        // N.B. the forbid coarsening will prevent me from being true if the block has been refined in the past
        status_lvl_ = (should_coarsen && status_lvl_ == M_ADAPT_SAME) ? M_ADAPT_COARSER : status_lvl_;
    // } else {
    //     m_log("we don't compute the details for this block, the computation is over: (%d && %d) || %d", forbid_coarsening, forbid_refinement, status_lvl_ == M_ADAPT_FINER);
    }

    // if (status_lvl_ == M_ADAPT_FINER) {
    //     m_log("block @ %f %f %f is to be refined", xyz_[0], xyz_[1], xyz_[2]);
    // }
    // if (status_lvl_ == M_ADAPT_COARSER) {
    //     m_log("block @ %f %f %f is to be coarsened", xyz_[0], xyz_[1], xyz_[2]);
    // }
    // prevent the blocks to have a none-determined status
    // status_lvl_ = (status_lvl_ == M_ADAPT_NONE) ? M_ADAPT_SAME : status_lvl_;
    //--------------------------------------------------------------------------
}

void GridBlock::UpdateStatusFromPatches(/* params */ const Wavelet* interp, std::list<Patch>* patch_list) {
    //--------------------------------------------------------------------------
    // m_assert(status_lvl_ == M_ADAPT_NONE, "trying to update a status which is already updated");
    //--------------------------------------------------------------------------
    // m_profStart(profiler, "patch");

    // prevent coarsening if we have finer neighbors
    // const bool forbid_coarsening = (local_children_.size() + ghost_children_.size()) > 0;
    // const bool forbid_refinement = (local_parent_.size() + ghost_parent_.size()) > 0;
    // const bool forbid_coarsening = ((local_children_.size() + ghost_children_.size()) > 0) || (level_ == 0);
    // const bool forbid_refinement = ((local_parent_.size() + ghost_parent_.size()) > 0) || (level_ == P8EST_QMAXLEVEL);

    // get the block length
    real_t len = p4est_QuadLen(this->level());

    // loop over the patches and determine if I am in it or not
    for (auto patch = patch_list->cbegin(); patch != patch_list->cend(); ++patch) {
        // if we already have the correct level or a higher one, we skip the patch
        if (this->level() > patch->level()) {
            // if not, we have a coarser block and we might want to refine if the location matches
            bool coarsen = true;// !forbid_coarsening;

            for (lda_t id = 0; id < 3; id++) {
                // we have to satisfy both the our max > min and the min < our max
                coarsen = coarsen &&
                          (this->xyz(id) < (patch->origin(id) + patch->length(id))) &&
                          (patch->origin(id) < (this->xyz(id) + len));
            }
            // register the status
            status_lvl_ = (coarsen) ? (M_ADAPT_COARSER) : status_lvl_;

            // if we found a matching patch, it's done
            if (coarsen) {
                // m_log("block @ %f %f %f coarsening! ", this->xyz(0), this->xyz(1), this->xyz(2));
                // m_profStop(profiler, "patch");
                return;
            }

        } else if (this->level() < patch->level()) {
            // if not, we have a coarser block and we might want to refine if the location matches
            bool refine = true;// !forbid_refinement;

            for (lda_t id = 0; id < 3; id++) {
                // we have to satisfy both the our max > min and the min < our max
                refine = refine &&
                         (this->xyz(id) < (patch->origin(id) + patch->length(id))) &&
                         (patch->origin(id) < (this->xyz(id) + len));
            }
            // register the status
            status_lvl_ = (refine) ? (M_ADAPT_FINER) : status_lvl_;

            // if we found a matching patch, it's done
            if (refine) {
                // m_log("block @ %f %f %f refining! ", this->xyz(0), this->xyz(1), this->xyz(2));
                // m_profStop(profiler, "patch");
                return;
            }
        }
    }
    // m_log("block @ %f %f %f not changing! ", this->xyz(0), this->xyz(1), this->xyz(2));
    // status_lvl_ = M_ADAPT_SAME;
    // m_assert(status_lvl_ != M_ADAPT_NONE, "the status of the block cannot be NONE");
    // m_profStop(profiler, "patch");
    return;
    //--------------------------------------------------------------------------
}

/**
 * @brief update the status to match the level requirements
 * 
 */
void GridBlock::UpdateStatusFromLevel(const level_t min_level, const level_t max_level) {
    //--------------------------------------------------------------------------
    // forbid to refine/coarsen if we are on the finest/coarsest level possible
    const bool forbid_coarsening = (level_ == min_level);
    const bool forbid_refinement = (level_ == max_level);

    status_lvl_ = (forbid_coarsening && status_lvl_ == M_ADAPT_COARSER) ? M_ADAPT_SAME : status_lvl_;
    status_lvl_ = (forbid_refinement && status_lvl_ == M_ADAPT_FINER) ? M_ADAPT_SAME : status_lvl_;
    //--------------------------------------------------------------------------
}

/**
 * @brief Forwards the refinement down the grid = if one of my finer neighbor wants to refine, I refine as well
 * 
 */
void GridBlock::UpdateStatusForwardRefinement() {
    //--------------------------------------------------------------------------
    // if one of my finer neighbor wants to refine, I have to refine as well
    bool force_refining = false;

    for (iblock_t icount = 0; icount < local_children_.size(); icount++) {
        force_refining = force_refining || (status_ngh_[M_LOC_CHILDREN][icount] == M_ADAPT_FINER);
    }
    for (iblock_t icount = 0; icount < ghost_children_.size(); icount++) {
        force_refining = force_refining || (status_ngh_[M_GLO_CHILDREN][icount] == M_ADAPT_FINER);
    }

    // if the block is SAME, we force the refinement
    // if the block is COARSER, we cannot coarsen because my neighbor needs me to refine
    // if the block is FINER, nothing changes
    status_lvl_ = (force_refining) ? M_ADAPT_FINER : status_lvl_;

    // if (force_refining) {
    //     m_log("block @ %f %f %f is forced to refine ",xyz_[0],xyz_[1],xyz_[2]);
    // }
    //--------------------------------------------------------------------------
}

/**
 * @brief update the status of the present block given its neighbor's status
 * 
 * for rule definitions, cfr the paper
 * 
 */
void GridBlock::UpdateStatusFromGlobalPolicy() {
    //--------------------------------------------------------------------------
    // forbid the coarsening if we have finer neighbors
    // forbid the refinement if I have coarser neighbors
    bool forbid_coarsening = (local_children_.size() + ghost_children_.size()) > 0;
    bool forbid_refining   = (local_parent_.size() + ghost_parent_.size()) > 0;

    // I cannot coarsen if one of my coarser or same level neighbor wants to refine -> rule (1)
    {
        for (iblock_t icount = 0; icount < local_parent_.size(); icount++) {
            forbid_coarsening = forbid_coarsening || (status_ngh_[M_LOC_PARENT][icount] == M_ADAPT_FINER);
        }
        for (iblock_t icount = 0; icount < ghost_parent_.size(); icount++) {
            forbid_coarsening = forbid_coarsening || (status_ngh_[M_GLO_PARENT][icount] == M_ADAPT_FINER);
        }
        for (iblock_t icount = 0; icount < local_sibling_.size(); icount++) {
            forbid_coarsening = forbid_coarsening || (status_ngh_[M_LOC_SIBLING][icount] == M_ADAPT_FINER);
        }
        for (iblock_t icount = 0; icount < ghost_sibling_.size(); icount++) {
            forbid_coarsening = forbid_coarsening || (status_ngh_[M_GLO_SIBLING][icount] == M_ADAPT_FINER);
        }
    }

    // if (forbid_refining && status_lvl_ == M_ADAPT_FINER) {
    //     m_log("block @ %f %f %f from %d to %d", xyz_[0], xyz_[1], xyz_[2], status_lvl_, M_ADAPT_SAME);
    // }

    // if (forbid_coarsening && status_lvl_ == M_ADAPT_COARSER) {
    //     m_log("block @ %f %f %f from %d to %d", xyz_[0], xyz_[1], xyz_[2], status_lvl_, M_ADAPT_SAME);
    // }

    // if I am forbidden from coarsening
    status_lvl_ = (forbid_coarsening && status_lvl_ == M_ADAPT_COARSER) ? M_ADAPT_SAME : status_lvl_;
    status_lvl_ = (forbid_refining && status_lvl_ == M_ADAPT_FINER) ? M_ADAPT_SAME : status_lvl_;

    //--------------------------------------------------------------------------
}

/**
 * @brief allocate the array for the different blocks needed
 * 
 * @param qid 
 * @param coarsen_vec 
 */
void GridBlock::SyncStatusInit() {
    //--------------------------------------------------------------------------

    // allocate the local status_ngh array
    iblock_t nblocks_total = (local_parent_.size() + ghost_parent_.size() +
                              local_sibling_.size() + ghost_sibling_.size() +
                              local_children_.size() + ghost_children_.size());

    // make sure of the order
    m_assert(M_LOC_PARENT == 0, "the first must be the Local parents");
    m_assert(M_GLO_PARENT == 1, "the first must be the Local parents");
    m_assert(M_LOC_SIBLING == 2, "the first must be the Local parents");
    m_assert(M_GLO_SIBLING == 3, "the first must be the Local parents");
    m_assert(M_LOC_CHILDREN == 4, "the first must be the Local parents");
    m_assert(M_GLO_CHILDREN == 5, "the first must be the Local parents");

    // split the array
    status_ngh_[M_LOC_PARENT]   = reinterpret_cast<short_t*>(m_calloc(nblocks_total * sizeof(short_t)));
    status_ngh_[M_GLO_PARENT]   = status_ngh_[M_LOC_PARENT] + local_parent_.size();
    status_ngh_[M_LOC_SIBLING] = status_ngh_[M_GLO_PARENT] + ghost_parent_.size();
    status_ngh_[M_GLO_SIBLING] = status_ngh_[M_LOC_SIBLING] + local_sibling_.size();
    status_ngh_[M_LOC_CHILDREN] = status_ngh_[M_GLO_SIBLING] + ghost_sibling_.size();
    status_ngh_[M_GLO_CHILDREN] = status_ngh_[M_LOC_CHILDREN] + local_children_.size();
    //--------------------------------------------------------------------------
}

void GridBlock::SyncStatusFill(const qid_t* qid, short_t* const coarsen_vec) {
    //--------------------------------------------------------------------------
    coarsen_vec[qid->cid] = (short_t)status_lvl_;
    //--------------------------------------------------------------------------
}

/**
 * @brief update my status based on the status from my finer neighbors
 * 
 * allocate the @ref status_siblings_neighbors_ array, which will be destroyed in the @ref SolveNeighbor function.
 */
void GridBlock::SyncStatusUpdate(const short_t* const status_vec, MPI_Win status_window) {
    //--------------------------------------------------------------------------
    const iblock_t n_coarser = local_parent_.size() + ghost_parent_.size();
    const iblock_t n_finer   = local_children_.size() + local_children_.size();
    const iblock_t n_local_coarser  = local_parent_.size();
    const iblock_t n_local_finer    = local_parent_.size();

    //..........................................................................
    // loop over the remote coarser blocks
    {
        iblock_t count = 0;
        for (auto* gblock : ghost_parent_) {
            // m_log("count = %d -> requestion block cum id = %ld at rank %d", count, displ, gblock->rank());
            m_assert(sizeof(short_t) == sizeof(short), "the two sizes must match to garantee mpi data types");
            MPI_Get(status_ngh_[M_GLO_PARENT] + count, 1, MPI_SHORT, gblock->rank(), gblock->cum_block_id(), 1, MPI_SHORT, status_window);
            ++count;
        }
    }
    {
        iblock_t count = 0;
        for (auto* gblock : ghost_children_) {
            // m_log("count = %d -> requestion block cum id = %ld at rank %d", count, displ, gblock->rank());
            m_assert(sizeof(short_t) == sizeof(short), "the two sizes must match to garantee mpi data types");
            MPI_Get(status_ngh_[M_GLO_CHILDREN] + count, 1, MPI_SHORT, gblock->rank(), gblock->cum_block_id(), 1, MPI_SHORT, status_window);
            ++count;
        }
    }
    {
        iblock_t count = 0;
        for (auto* gblock : ghost_sibling_) {
            // m_log("count = %d -> requestion block cum id = %ld at rank %d", count, displ, gblock->rank());
            m_assert(sizeof(short_t) == sizeof(short), "the two sizes must match to garantee mpi data types");
            MPI_Get(status_ngh_[M_GLO_SIBLING] + count, 1, MPI_SHORT, gblock->rank(), gblock->cum_block_id(), 1, MPI_SHORT, status_window);
            ++count;
        }
    }
    //..........................................................................
    // get the local ones now
    {
        iblock_t count = 0;
        for (auto* gblock : local_parent_) {
            status_ngh_[M_LOC_PARENT][count] = status_vec[gblock->cum_block_id()];
            ++count;
        }
    }

    {
        iblock_t count = 0;
        for (auto* gblock : local_children_) {
            status_ngh_[M_LOC_CHILDREN][count] = status_vec[gblock->cum_block_id()];
            ++count;
        }
    }
    {
        iblock_t count = 0;
        for (auto* gblock : local_sibling_) {
            status_ngh_[M_LOC_SIBLING][count] = status_vec[gblock->cum_block_id()];
            ++count;
        }
    }
    //--------------------------------------------------------------------------
}

void GridBlock::SyncStatusFinalize() {
    //---------------------------------------------------------------------------
    m_free(status_ngh_[M_LOC_PARENT]);
    //---------------------------------------------------------------------------
}
void GridBlock::SmoothResolutionJump(const Wavelet* interp, std::map<std::string, Field*>::const_iterator field_start, std::map<std::string, Field*>::const_iterator field_end){
    // the status level has to be 0, otherwise it means that one of the block is not coarsened
    // m_assert(status_lvl_ != M_ADAPT_NONE, "we should have made a decision here");
    //--------------------------------------------------------------------------
    // m_profStart(profiler, "smooth jump");
    // reset the temp memory to 0.0
    // memset(coarse_ptr_(), 0, CartBlockMemNum(1) * sizeof(real_t));
    // data_ptr mask      = coarse_ptr_(0, this);
    // real_t*  mask_data = mask.Write();
    // reset to 0
    coarse_ptr_.MemSetZero();
    MemLayout myself = BlockLayout();
    MemData   mask(&coarse_ptr_,&myself);

    // collect the sizes as "by-copy capture of value of abstract type 'const Wavelet' is not allowed"
    const bidx_t n_criterion_front = interp->ndetail_citerion_extend_front();
    const bidx_t n_criterion_back = interp->ndetail_citerion_extend_back();
    const bidx_t n_smooth_front = interp->ndetail_smooth_extend_front();
    const bidx_t n_smooth_back =interp->ndetail_smooth_extend_back();

    //................................................
    // lambda to obtain the smoothing pattern
    auto mask_smooth = [=](const short_t status_ngh, const iface_t ibidule) -> void {
        // create the lambda to put 1.0
        auto set_mask_to_one = [=](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
            mask(i0, i1, i2) = 1.0;
        };
        // if the neighbor is a newly created block -> smooth
        if (status_ngh == M_ADAPT_NEW_COARSE) {
            m_assert(this->status_level() <= M_ADAPT_SAME, "if my coarser neighbor has been newly created, I cannot have something different than SAME (now %d)", status_lvl_);

            // get the sign of the ibidule
            real_t sign[3];
            GhostGetSign(ibidule, sign);

            // create the start and end indexes
            MemSpan smooth_span;
            MemSpan block_span = this->BlockSpan();
            
#pragma unroll 3
            for (lda_t ida = 0; ida < 3; ++ida) {
                if (sign[ida] > 0.5) {
                    // // my ngh assumed 0 details in my block
                    // smooth_start[ida] = this->end(ida) - interp->ndetail_citerion_extend_front();
                    // // the number of my ngh details influencing my values
                    // smooth_end[ida] = this->end(ida) + interp->ndetail_smooth_extend_back();
                    // my ngh assumed 0 details in my block
                    smooth_span.start[ida] = block_span.end[ida] - n_criterion_front;
                    smooth_span.end[ida]   = block_span.end[ida] + n_smooth_back;
                } else if (sign[ida] < (-0.5)) {
                    // // my ngh assumed 0 details in my block
                    // smooth_start[ida] = this->start(ida) - interp->ndetail_smooth_extend_front();
                    // // the number of my ngh details influencing my values
                    // smooth_end[ida] = this->start(ida) + interp->ndetail_citerion_extend_back();
                    // my ngh assumed 0 details in my block
                    smooth_span.start[ida] = block_span.start[ida] - n_smooth_front;
                    smooth_span.end[ida]   = block_span.start[ida] + n_criterion_back;
                } else {
                    // even in the directions orthogonal to ibidule, the details must be killed!
                    // as my neighbor, which might be fine will kill them as well
                    // smooth_start[ida] = this->start(ida) - interp->ndetail_smooth_extend_front();
                    // smooth_end[ida]   = this->end(ida) + interp->ndetail_smooth_extend_back();
                    smooth_span.start[ida] = block_span.start[ida] - n_smooth_front;
                    smooth_span.end[ida]   = block_span.end[ida] + n_smooth_back;
                }
            }
            // m_log("ibidule = %d -> maks = 1.0 from %d %d %d to %d %d %d", ibidule,
            //       smooth_start[0], smooth_start[1], smooth_start[2],
            //       smooth_end[0], smooth_end[1], smooth_end[2]);
            // apply it
            // for_loop(&set_mask_to_one, smooth_start, smooth_end);
            for_loop(&set_mask_to_one,&smooth_span);
        }
    };

    // for each ghost block, set the mask to 1.0 if needed
    iblock_t block_count = 0;
    for (auto* gblock : local_parent_) {
        // if ((status_ngh_[block_count] != M_ADAPT_SAME) && (this->status_level() != M_ADAPT_SAME)) {
        //     m_log("block thinks that his neigbor %d (count=%d) has been modified while the block has been changed", gblock->cum_block_id(), block_count);
        //     m_assert(false, "oouuups");
        // }
        mask_smooth(status_ngh_[M_LOC_PARENT][block_count], gblock->ibidule());
        // update the counter
        ++block_count;
    }
    m_assert(block_count == local_parent_.size(), "the two numbers must match: %d vs %ld", block_count, local_parent_.size());
    
    block_count = 0;
    for (auto* gblock : ghost_parent_) {
        // if ((status_ngh_[block_count] != M_ADAPT_SAME) && (this->status_level() != M_ADAPT_SAME)) {
        //     m_log("block thinks that his neigbor %d (count=%d) has been modified while the block has been changed", gblock->cum_block_id(), block_count);
        //     m_assert(false, "oouuups");
        // }
        mask_smooth(status_ngh_[M_GLO_PARENT][block_count], gblock->ibidule());
        // update the counter
        ++block_count;
    }
    m_assert(block_count == ghost_parent_.size(), "the two numbers must match: %d vs %ld", block_count, ghost_parent_.size());

    //................................................
    // smooth depending on the mask
    // SubBlock block_src(this->gs(), this->stride(), -interp->nghost_front(), M_N + interp->nghost_back());
    // SubBlock block_det(this->gs(), this->stride(), -interp->ndetail_smooth_extend_front(), M_N + interp->ndetail_smooth_extend_back());
    MemSpan block_span = this->BlockSpan();
    MemSpan span_src(-interp->nghost_front(), M_N + interp->nghost_back());
    MemSpan span_det(-interp->ndetail_smooth_extend_front(), M_N + interp->ndetail_smooth_extend_back());

    // do it for every field
    for (auto fid = field_start; fid != field_end; ++fid) {
        auto* current_field = fid->second;
        if (!current_field->is_temp()) {
            for (lda_t ida = 0; ida < current_field->lda(); ida++) {
                // interp->SmoothOnMask(&block_src, this, this->data(current_field, ida), &block_det, mask);
                const MemData data = this->data(current_field,ida);
                interp->SmoothOnMask(&span_src,&block_span,&data,&span_det,&mask);//&block_src, this, this->data(current_field, ida), &block_det, mask);
            }
        }
    }

    // free the status array
    // m_free(status_ngh_);
    // m_profStop(profiler, "smooth jump");
    //--------------------------------------------------------------------------
}

/**
 * @brief computes the Maximum and Minimun detail coefficient on the block
 * 
 * the function updates the ongoing maxmin array if needed (no initialization is performed in the function)
 * 
 * @param qid the quadrant ID of the block
 * @param interp the Wavelet
 * @param criterion the criterion field, must be scalar (for now)
 * @param maxmin the minmax that is used to collect on the blocks
 * @param max_blocks an array containing the max of every block, can be nullptr
 */
void GridBlock::MaxMinDetails(const Wavelet* interp, const Field* criterion, real_t maxmin[2],
                              bidx_t* max_blocks, const real_t max_cat, const real_t min_cat, const short_t n_cat) {
    m_assert(criterion->lda() == 1, "field must be a scalar");
    //-------------------------------------------------------------------------
    real_t block_maxmin[2] = {0.0, 0.0};

    const MemSpan span_src(-interp->nghost_front(), M_N + interp->nghost_back());
    const MemSpan span_det = this->BlockSpan();
    ConstMemData  data_src = this->ConstData(criterion, 0);
    MemData       data_det(nullptr);
    interp->Details(&span_src, &data_src, &span_det, &data_det, 0.0, block_maxmin);

    maxmin[0] = m_max(maxmin[0], block_maxmin[0]);
    maxmin[1] = m_min(maxmin[1], block_maxmin[1]);

    // we store the categories
    if (max_blocks != nullptr) {
        // get the category split
        const real_t h_cat  = (log10(max_cat) - log10(min_cat)) / n_cat;
        short_t      id_cat = m_max(0, m_min(n_cat - 1, (log10(block_maxmin[0]) - log10(min_cat)) / h_cat));

        // got
        m_assert(n_cat > id_cat && id_cat >= 0, "the cat id = %d must be >=0 and < %d", id_cat, n_cat);
        max_blocks[id_cat] += 1;
    }
    //--------------------------------------------------------------------------
}

/**
 * @brief Compute the detail coefficient and store them in the field details
 * 
 * @param interp the wavelet object
 * @param criterion the criterion field
 * @param details the detail field with the compute detail values
 */
void GridBlock::StoreDetails(const Wavelet* interp, const Field* criterion, const Field* details) {
    const bidx_t ghost_len_interp[2] = {interp->nghost_front(), interp->nghost_back()};
    m_assert(criterion->ghost_status(ghost_len_interp), "the field <%s> must have up-to-date ghosts", criterion->name().c_str());
    m_assert(criterion->lda() == details->lda(), "field <%s> and <%s> must have the same size", criterion->name().c_str(), details->name().c_str());
    //--------------------------------------------------------------------------
    for (lda_t ida = 0; ida < criterion->lda(); ida++) {
        // SubBlock block_src(this->gs(), this->stride(), -interp->nghost_front(), M_N + interp->nghost_back());
        // interp->StoreDetails(&block_src, this->data(criterion, ida), this, this->data(details, ida));
        const MemSpan      span_src(-interp->nghost_front(), M_N + interp->nghost_back());
        const MemSpan      me       = BlockSpan();
        const ConstMemData data_src = this->ConstData(criterion, ida);
        const MemData      data_trg = this->data(details, ida);
        interp->StoreDetails(&span_src, &data_src, &me, &data_trg);
    }
    //--------------------------------------------------------------------------
}


/**
 * @brief resolve the dependency list created while adapting the mesh (see cback_UpdateDependency() ) by interpolating the needed blocks
 * 
 * @param interp the Wavelet object to use for the interpolation/refinement
 * @param field_start the first field to take into account
 * @param field_end the last field to take into account
 */
void GridBlock::SolveDependency(const Wavelet* interp, std::map<std::string, Field*>::const_iterator field_start, std::map<std::string, Field*>::const_iterator field_end) {
    m_assert(n_dependency_active_ == 0 || n_dependency_active_ == 1 || n_dependency_active_ == P8EST_CHILDREN, "wrong value for n_dependency_active_");
    // m_profStart(profiler, "solve dependency");
    //--------------------------------------------------------------------------
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
        MemSpan     span_src(src_start, src_end);
        // for every field on my parent, interpolate it
        m_assert(mem_map_.empty(), "the block should be empty here");
        for (auto fid = field_start; fid != field_end; ++fid) {
            auto current_field = fid->second;
            // allocate the field on me (has not been done yet)
            if (!current_field->is_expr()) {
                this->AddField(current_field);
                // refine for every dimension, if the field is not temp
                if (!current_field->is_temp()) {
                    // check the ghost ability
                    const bidx_t ghost_len[2] = {interp->nghost_front_refine(), interp->nghost_back_refine()};
                    m_assert(current_field->ghost_status(ghost_len), "The field <%s> must have enough valid GP for the refinement - required %d %d, known %d %d", current_field->name().c_str(), ghost_len[0], ghost_len[1], current_field->get_ghost_len(0), current_field->get_ghost_len(1));
                    for (sid_t ida = 0; ida < current_field->lda(); ida++) {
                        // get the pointers
                        // interp->Interpolate(-1, shift, &mem_src, root->data(current_field, ida), this, this->data(current_field, ida));
                        // MemData
                        const ConstMemData data_src = root->ConstData(current_field, ida);
                        const MemSpan      span_trg = this->BlockSpan();
                        const MemData      data_trg = this->data(current_field, ida);
                        interp->Interpolate(-1, shift, &span_src, &data_src, &span_trg, &data_trg);
                    }
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
        m_assert(this->status_level() == M_ADAPT_NEW_COARSE || this->status_level() == M_ADAPT_SAME, "my status must be M_ADAPT_COARSER or M_ADAPT_SAME instead of %d", this->status_level());
        // I have 8 deps, I am a root, waiting data from coarsening of my children
        //allocate the new fields
        m_assert(mem_map_.size() == 0, "the block should be empty here");
        for (auto fid = field_start; fid != field_end; ++fid) {
            auto current_field = fid->second;
            // allocate the field on me (has not been done yet)
            if (!current_field->is_expr()) {
                this->AddField(current_field);
            }
        }

        // I am a parent and I need to fillout my children
        for (sid_t childid = 0; childid < P8EST_CHILDREN; ++childid) {
            GridBlock* child_block = this->PopDependency(childid);
            m_assert(child_block->level() - this->level() == 1, "the child block is not a child");
            m_assert(childid == p4est_GetChildID(child_block->xyz(), child_block->level()), "the two ids must match");
            m_assert(child_block->status_level() == M_ADAPT_COARSER, "the status of the new root must be M_ADAPT_NEW_COARSE instead of %d", child_block->status_level());

            // get the shift for me = child
            const lid_t shift[3]     = {-M_N * ((childid % 2)), -M_N * ((childid % 4) / 2), -M_N * ((childid / 4))};
            const lid_t trg_start[3] = {M_NCENTER * ((childid % 2)), M_NCENTER * ((childid % 4) / 2), M_NCENTER * ((childid / 4))};
            const lid_t trg_end[3]   = {trg_start[0] + M_NCENTER, trg_start[1] + M_NCENTER, trg_start[2] + M_NCENTER};
            MemSpan     span_trg(trg_start, trg_end);
            // SubBlock    mem_trg(M_GS, M_STRIDE, trg_start, trg_end);
            // and an extended source block for my child
            // const lid_t src_start[3] = {-M_GS, -M_GS, -M_GS};
            // const lid_t src_end[3]   = {M_N + M_GS, M_N + M_GS, M_N + M_GS};
            // MemSpan     span_src(src_start, src_end);
            // SubBlock    mem_src(M_GS, M_STRIDE, src_start, src_end);

            // for every field, we interpolate it
            for (auto fid = field_start; fid != field_end; ++fid) {
                auto current_field = fid->second;
                if (!current_field->is_temp()) {
                    const bidx_t ghost_len[2] = {interp->nghost_front_coarsen(), interp->nghost_back_coarsen()};
                    m_assert(current_field->ghost_status(ghost_len), "The field must have enough valid GP for the refinement - required %d %d, known %d %d", ghost_len[0], ghost_len[1], current_field->get_ghost_len(0), current_field->get_ghost_len(1));
                    MemSpan      span_src = child_block->ExtendedSpan(ghost_len);
                    // interpolate for every dimension
                    for (sid_t ida = 0; ida < current_field->lda(); ida++) {
                        // get the pointers
                        const ConstMemData data_child = child_block->ConstData(current_field, ida);
                        const MemData      data_root  = this->data(current_field, ida);
                        interp->Interpolate(1, shift, &span_src, &data_child, &span_trg, &data_root);
                        // interp->Interpolate(1, shift, &mem_src, child_block->data(current_field, ida), &mem_trg, this->data(current_field, ida));
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
    // m_profStop(profiler, "solve dependency");
    //--------------------------------------------------------------------------
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
    //--------------------------------------------------------------------------
    --n_dependency_active_;
    GridBlock* block      = dependency_[child_id];
    dependency_[child_id] = nullptr;
    return block;
    //--------------------------------------------------------------------------
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
    //--------------------------------------------------------------------------
    ++n_dependency_active_;
    dependency_[child_id] = dependent_block;
    //--------------------------------------------------------------------------
}

/**
 * @brief build the list of ghosts
 * 
 * @param qid the current quadrant ID
 * @param grid the grid to use to recover the ghosts etc
 * @param interp the wavelet to use (for coarse ghost sizes etc)
 * @param local2disp_window the displacement information for RMA
 */
void GridBlock::GhostInitLists(const qid_t* qid, const ForestGrid* grid, const Wavelet* interp, MPI_Win local2disp_window) {
    //--------------------------------------------------------------------------
    // allocate the ghost pointer, which is reused for the wavelets smoothing
    // size_t alloc_size = m_max(interp->CoarseSize(), m_blockmemsize(1));
    // AllocateCoarsePtr(alloc_size);
    // m_log("I allocate %ld doubles",alloc_size);
    // m_assert(interp->CoarseSize() <= CartBlockMemNum(1), "the coarse size must be smaller than a blockmemsize to fit in the coarse memory");
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
    // the ghost points are computed based on the maximal number of ghosts possible!
    // lid_t  block_min[3], block_max[3];
    real_t block_len[3];
    real_t coarse_hgrid[3];
    for (lda_t id = 0; id < 3; id++) {
        m_assert(level() >= 0, "the level=%d must be >=0", level());
        // set the number of ghost to compute
        // block_min[id]    = -interp->nghost_front();
        // block_max[id]    = M_N + interp->nghost_back();
        block_len[id]    = p4est_QuadLen(level());
        coarse_hgrid[id] = CoarseHGrid(p4est_QuadLen(level()));
    }
    const bidx_t block_ghost_len[2] = {M_GS, M_GS};
    const bidx_t block_core_len     = M_N;

    rank_t my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    // get the current gs and stride
    // const bidx_t block_gs     = gs();
    // const bidx_t block_stride = stride();
    MemLayout block_layout = BlockLayout();

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
                const MemSpan me = BlockSpan();
                PhysBlock* pb = new PhysBlock(&ghost_len_, ibidule, &me);  //, interp->nghost_front(), interp->nghost_back());
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
            m_verb("reading the list: adress: %p  and rank %d -> is ghost? %d", nghq, ngh_rank, isghost);

            // get the sign, i.e. the normal to the face, the edge of the corner we consider
            real_t sign[3];
            GhostGetSign(ibidule, sign);

            //................................................
            // get the position of the neighbor, as seen by me!!! may be different than the actual one if there is a periodic bc
            real_t ngh_pos[3];
            if (!isghost) {
                // cannot use the p8est function because the which_tree is not assigned, so we retrieve the position through the block
                // GridBlock* ngh_block = *(reinterpret_cast<GridBlock**>(nghq->p.user_data));
                GridBlock* ngh_block = p4est_GetGridBlock(nghq);
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
                    // GBLocal* gb = new GBLocal(block_gs, block_stride, &ghost_len_, ibidule, ngh_cum_id);
                    GBLocal* gb = new GBLocal(&ghost_len_, ibidule, ngh_cum_id);
                    gb->Intersect(/* source */ ngh_block->level(), ngh_pos, ngh_hgrid, ngh_len,
                                  /* traget */ level_, xyz_, hgrid_, block_ghost_len, block_core_len);  // block_min, block_max);
                    gb->data_src(ngh_block);

                    //#pragma omp critical
                    local_sibling_.push_back(gb);

                } else if (nghq->level < level()) {
                    // m_log("creating a coarser");
                    // parent: source = neighbor, target = me
                    GBLocal* gb = new GBLocal(&ghost_len_, ibidule, ngh_cum_id);
                    gb->Intersect(/* source */ ngh_block->level(), ngh_pos, ngh_hgrid, ngh_len,
                                  /* target */ level_, xyz_, hgrid_, block_ghost_len, block_core_len);  // block_min, block_max);
                    gb->data_src(ngh_block);
                    //#pragma omp critical
                    local_parent_.push_back(gb);

                    // m_log("creating a reverse");
                    // the children: the source = the coarse myself, target = my neighbor
                    GBLocal* invert_gb = new GBLocal( &ghost_len_, -1, ngh_cum_id);
                    invert_gb->Intersect(/* source */ level_ - 1, xyz_, coarse_hgrid, block_len,
                                         /* target */ ngh_block->level(), ngh_pos, ngh_hgrid, block_ghost_len, block_core_len);  // block_min, block_max);
                    invert_gb->data_src(ngh_block);
                    //#pragma omp critical
                    local_parent_reverse_.push_back(invert_gb);
                } else if (nghq->level > level()) {
                    // m_log("creating a finer");
                    m_assert((nghq->level - level_) == 1, "The delta level is not correct: %d - %d", nghq->level, level_);
                    // register the coarse
                    GBLocal* gb = new GBLocal(&ghost_len_, ibidule, ngh_cum_id);
                    //#pragma omp critical
                    local_children_.push_back(gb);
                } else {
                    m_assert(false, "The delta level is not correct: %d - %d", nghq->level, level());
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
                    GBMirror* gb = new GBMirror(&ghost_len_, ibidule, ngh_cum_id, ngh_rank);
                    gb->Intersect(/* source */ nghq->level, ngh_pos, ngh_hgrid, ngh_len,
                                  /* target */ level(), xyz(), hgrid(), block_ghost_len, block_core_len);  // block_min, block_max);
                    // ask the displacement (will be available later, when completing the call)
                    m_assert(sizeof(bidx_t) == sizeof(int), "the two sizes must match");
                    MPI_Get(gb->data_src_ptr(), 1, MPI_AINT, ngh_rank, ngh_cum_id, 1, MPI_AINT, local2disp_window);
                    //#pragma omp critical
                    ghost_sibling_.push_back(gb);

                }
                //................................................
                else if (nghq->level < level_) {
                    // parent: source = neighbor, target = me
                    GBMirror* gb = new GBMirror(&ghost_len_, ibidule, ngh_cum_id, ngh_rank);
                    gb->Intersect(/* source */ nghq->level, ngh_pos, ngh_hgrid, ngh_len,
                                  /* target */ level(), xyz(), hgrid(), block_ghost_len, block_core_len);  // block_min, block_max);
                    // ask the displacement (will be available later, when completing the call
                    MPI_Get(gb->data_src_ptr(), 1, MPI_AINT, ngh_rank, ngh_cum_id, 1, MPI_AINT, local2disp_window);
                    //#pragma omp critical
                    ghost_parent_.push_back(gb);

                    // I compute my own contribution to my neighbor ghost points
                    GBMirror* invert_gb = new GBMirror(&ghost_len_, ibidule, ngh_cum_id, ngh_rank);
                    invert_gb->Intersect(/* source */ level() - 1, xyz(), coarse_hgrid, block_len,
                                         /* target */ nghq->level, ngh_pos, ngh_hgrid, block_ghost_len, block_core_len);  // block_min, block_max);
                    MPI_Get(invert_gb->data_src_ptr(), 1, MPI_AINT, ngh_rank, ngh_cum_id, 1, MPI_AINT, local2disp_window);
                    //#pragma omp critical
                    ghost_parent_reverse_.push_back(invert_gb);
                }
                //................................................
                else if (nghq->level > level_) {
                    const real_t ngh_hgrid_coarse[3] = {CoarseHGrid(ngh_len[0]), CoarseHGrid(ngh_len[1]), CoarseHGrid(ngh_len[2])};

                    // children: source = coarse version of my neighbor, target = myself
                    GBMirror* gb = new GBMirror(&ghost_len_, ibidule, ngh_cum_id, ngh_rank);
                    gb->Intersect(/* source */ nghq->level - 1, ngh_pos, ngh_hgrid_coarse, ngh_len,
                                  /* target */ level(), xyz(), hgrid(), block_ghost_len, block_core_len);  // block_min, block_max);
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
    // //--------------------------------------------------------------------------
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
    //--------------------------------------------------------------------------
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
    //--------------------------------------------------------------------------
}

/**
 * @brief update my ghosts lists to have the exact number of ghost required by the user
 * 
 * @param ghost_sizes the ghost size in front and at the back of the block
 */
void GridBlock::GhostUpdateSize(const bidx_t ghost_len[2]) {
    m_assert(ghost_len[0] >= 0, "the ghost size = %d must be >=0", ghost_len[0]);
    m_assert(ghost_len[1] >= 0, "the ghost size = %d must be >=0", ghost_len[1]);
    //--------------------------------------------------------------------------
    // store them for me
    ghost_len_[0] = ghost_len[0];
    ghost_len_[1] = ghost_len[1];

    // // get the lambda to execute on each list
    // auto adapt_len = [ghost_len](auto& list) -> void {
    //     for (GhostBlock* __restrict block : list) {
    //         block->ghost_len(ghost_len);
    //     }
    // };

    // // local
    // adapt_len(local_sibling_);
    // adapt_len(local_parent_);
    // adapt_len(local_children_);
    // adapt_len(local_parent_reverse_);
    // // ghost
    // adapt_len(ghost_sibling_);
    // adapt_len(ghost_parent_);
    // adapt_len(ghost_children_);
    // adapt_len(ghost_parent_reverse_);
    // // physics
    // adapt_len(phys_);
    //--------------------------------------------------------------------------
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
void GridBlock::GhostGet_Cmpt(const Field* field, const lda_t ida, const Wavelet* interp) {
    //--------------------------------------------------------------------------
    // get the siblings
    {
        // const SubBlock bsrc_neighbor(M_GS, M_STRIDE, 0, M_N);
        // m_assert(bsrc_neighbor.core() == M_N, "the core = %d must be %d", bsrc_neighbor.core(), M_N);
        // get the copy of local siblings
        for (auto* const gblock : local_sibling_) {
            const CartBlock*   ngh_block = gblock->data_src();
            const MemSpan      span_src  = gblock->SourceSpan();
            const MemSpan      span_trg  = gblock->GetSpan();
            const MemData      data_trg  = this->data(field, ida);
            const ConstMemData data_src  = ngh_block->ConstData(field, ida);
            // copy the information
            interp->Copy(gblock->dlvl(), gblock->shift(), &span_src, & data_src, & span_trg, &data_trg);
            // interp->Copy(gblock->dlvl(), gblock->shift(), &bsrc_neighbor, data_src, block_trg, data_trg);
        }
    }

    //................................................
    const bool do_coarse = (local_parent_.size() + ghost_parent_.size()) > 0;
    if (do_coarse) {
        // we use the full neighbor span
        // const SubBlock bsrc_neighbor(M_GS, M_STRIDE, 0, M_N);
        // m_assert(bsrc_neighbor.core() == M_N, "the core = %d must be %d", bsrc_neighbor.core(), M_N);
        const MemLayout block_layout  = BlockLayout();
        const MemLayout coarse_layout = CoarseLayout(interp);

        auto copy_2_coarse = [=,&interp](const GBLocal* gblock) -> void {
            // get the neighbor info = source
            const CartBlock*   ngh_block = gblock->data_src();
            const MemSpan      span_ngh  = gblock->SourceSpan();
            const ConstMemData data_ngh  = ngh_block->ConstData(field, ida);

            // get the coarse span = target
            MemSpan span_coarse;
            gblock->GetCoarseSpan(&block_layout, &coarse_layout, &span_coarse);
            const MemData data_coarse(&coarse_ptr_,&coarse_layout);

            interp->Copy(gblock->dlvl() + 1, gblock->shift(), &span_ngh, &data_ngh, &span_coarse,& data_coarse);
        };

        for (auto* const gblock : local_sibling_){
            m_assert((gblock->dlvl() + 1) == 1, "the difference of level MUST be 1 --> dlvl == %d", (gblock->dlvl() + 1));
            copy_2_coarse(gblock);
        }
        for (auto* const gblock : local_parent_) {
            m_assert((gblock->dlvl() + 1) == 0, "the difference of level MUST be 0 --> dlvl == %d", (gblock->dlvl() + 1));
            copy_2_coarse(gblock);
        }

        // // copy the siblings to coarse
        // for (auto* const gblock : local_sibling_) {
        //     // get the neighbor info = source
        //     const GridBlock* ngh_block = gblock->data_src();
        //     const MemSpan span_ngh = ngh_block->BlockSpan();
        //     const ConstMemData data_ngh = ngh_block->data(field, ida);

        //     // get the coarse span = target
        //     MemSpan span_coarse;
        //     gblock->GetCoarseSpan(block_layout,coarse_layout,&span_coarse);
        //     const MemData      data_coarse(coarse_ptr_, coarse_layout);
            
        //     // SubBlock   block_trg;
        //     // // interp->CoarseFromFine(gblock, &block_trg);
        //     // gblock->ToCoarse(interp, &block_trg);
        //     // data_ptr data_src = ngh_block->data(field, ida);
        //     // data_ptr data_trg = coarse_ptr_(0, &block_trg);
        //     // interpolate, the level is 1 coarser and the shift is unchanged
        //     m_assert((gblock->dlvl() + 1) == 1, "the difference of level MUST be 1");
        //     interp->Copy(gblock->dlvl()+1,gblock->shift(),span_ngh,data_ngh,span_coarse,data_coarse);
        //     // interp->Copy(gblock->dlvl() + 1, gblock->shift(), &bsrc_neighbor, data_src, &block_trg, data_trg);
        // }
        // // copy the parents to coarse
        // for (auto* const gblock : local_parent_) {
        //     GridBlock* ngh_block = gblock->data_src();
        //     SubBlock   block_trg;
        //     // interp->CoarseFromFine(gblock, &block_trg);
        //     gblock->ToCoarse(interp, &block_trg);
        //     data_ptr data_src = ngh_block->data(field, ida);
        //     data_ptr data_trg = coarse_ptr_(0, &block_trg);
        //     // interpolate, the level is 1 coarser and the shift is unchanged
        //     m_assert((gblock->dlvl() + 1) == 0, "the difference of level MUST be 1");
        //     interp->Copy(gblock->dlvl() + 1, gblock->shift(), &bsrc_neighbor, data_src, &block_trg, data_trg);
        // }
    }
    //--------------------------------------------------------------------------
}

/**
 * @brief Do the first part of GhostGet, everything that requires MPI
 * 
 * @param field 
 * @param ida 
 * @param interp 
 * @param mirrors_window 
 */
void GridBlock::GhostGet_Post(const Field* field, const lda_t ida, const Wavelet* interp, MPI_Win mirrors_window) {
    //--------------------------------------------------------------------------
    // get the siblings
    {
        // const SubBlock bsrc_neighbor(M_GS, M_STRIDE, 0, M_N);
        // const data_ptr data_trg = data(field, ida);

        // // start the RMA for the sibling neighbors
        for (auto* const gblock : ghost_sibling_) {
            // source
            const MPI_Aint  disp_src   = gblock->data_src();
            const rank_t    rank_src   = gblock->rank();
            const MemLayout layout_src = gblock->SourceLayout();
            const MemSpan   span_src   = gblock->SourceSpan();
            // target
            const MemSpan   span_trg   = gblock->GetSpan();
            const MemLayout layout_trg = this->BlockLayout();
            const MemData   data_trg   = this->data(field, ida);
            // copy the information
            GetRma(gblock->dlvl(), gblock->shift(), &layout_src, &span_src, disp_src, rank_src, &layout_trg, &span_trg, &data_trg, mirrors_window);
            //     const MPI_Aint   disp_src  = gblock->data_src();
            //     const rank_t     disp_rank = gblock->rank();
            //     const MemSpan* block_trg = gblock;
            //     // copy the information
            //     interp->GetRma(gblock->dlvl(), gblock->shift(), &bsrc_neighbor, disp_src, block_trg, data_trg, disp_rank, mirrors_window);
        }
    }
    //................................................
    const bool do_coarse = (local_parent_.size() + ghost_parent_.size()) > 0;
    if (do_coarse) {
        coarse_ptr_.MemSetZero();
        MemLayout layout_coarse = this->CoarseLayout(interp);

        auto rma_2_coarse = [=](const GBMirror* gblock) {
            // source = block on the other side
            MemLayout layout_src = gblock->SourceLayout();
            MemSpan   span_src   = gblock->SourceSpan();
            MPI_Aint  disp_src   = gblock->data_src();
            rank_t    rank_src   = gblock->rank();

            // target = coarse pointer
            MemSpan         span_coarse;
            const MemLayout layout = this->BlockLayout();
            gblock->GetCoarseSpan(&layout, &layout_coarse, &span_coarse);
            MemData data_trg(&coarse_ptr_, &layout_coarse);
            GetRma((gblock->dlvl() + 1), gblock->shift(), &layout_src, &span_src, disp_src, rank_src, &layout_coarse, &span_coarse, &data_trg, mirrors_window);
        };

        for (auto* const gblock : ghost_sibling_) {
            m_assert((gblock->dlvl() + 1) == 1, "the difference of level MUST be 1 --> dlvl == %d", (gblock->dlvl() + 1));
            rma_2_coarse(gblock);
        }

        for (auto* const gblock : ghost_parent_) {
            m_assert((gblock->dlvl() + 1) == 0, "the difference of level MUST be 0 --> dlvl == %d", (gblock->dlvl() + 1));
            rma_2_coarse(gblock);
        }

        //         // interp->GetRma((gblock->dlvl() + 1), gblock->shift(), &bsrc_neighbor, disp_src, &block_trg, data_trg, disp_rank, mirrors_window);

        //         // // get the associated coarse block
        //         // SubBlock block_trg;
        //         // // interp->CoarseFromFine(gblock, &block_trg);
        //         // gblock->ToCoarse(interp, &block_trg);
        //         // data_ptr data_trg = coarse_ptr_(0, &block_trg);  // + m_zeroidx(0, &block_trg);
        //         // // interpolate, the level is 1 coarser and the shift is unchanged
        //         // m_assert((gblock->dlvl() + 1) == 1, "the difference of level MUST be 1 or 0");
        //         // interp->GetRma((gblock->dlvl() + 1), gblock->shift(), &bsrc_neighbor, disp_src, &block_trg, data_trg, disp_rank, mirrors_window);

        //     // reset the coarse memory
        //     // memset(coarse_ptr_(), 0, interp->CoarseSize(ghost_len_) * sizeof(real_t));

        //     // we use the full neighbor span
        //     // const SubBlock bsrc_neighbor(M_GS, M_STRIDE, 0, M_N);
        //     //................................................
        //     // RMA the sibligns ghosts to the tmp
        //     for (auto* const gblock : ghost_sibling_) {
        //         MPI_Aint disp_src  = gblock->data_src();
        //         rank_t   disp_rank = gblock->rank();
        //         // get the associated coarse block
        //         SubBlock block_trg;
        //         // interp->CoarseFromFine(gblock, &block_trg);
        //         gblock->ToCoarse(interp, &block_trg);
        //         data_ptr data_trg = coarse_ptr_(0, &block_trg);  // + m_zeroidx(0, &block_trg);
        //         // interpolate, the level is 1 coarser and the shift is unchanged
        //         m_assert((gblock->dlvl() + 1) == 1, "the difference of level MUST be 1 or 0");
        //         interp->GetRma((gblock->dlvl() + 1), gblock->shift(), &bsrc_neighbor, disp_src, &block_trg, data_trg, disp_rank, mirrors_window);
        //     }
        //     // RMA the parent ghosts to the tmp
        //     for (auto* const gblock : ghost_parent_) {
        //         MPI_Aint disp_src  = gblock->data_src();
        //         rank_t   disp_rank = gblock->rank();
        //         // get the associated coarse block
        //         SubBlock block_trg;
        //         // interp->CoarseFromFine(gblock, &block_trg);
        //         gblock->ToCoarse(interp, &block_trg);
        //         data_ptr data_trg = coarse_ptr_(0, &block_trg);
        //         // interpolate, the level is 1 coarser and the shift is unchanged
        //         m_assert((gblock->dlvl() + 1) == 0, "the difference of level MUST be 1 or 0");
        //         interp->GetRma((gblock->dlvl() + 1), gblock->shift(), &bsrc_neighbor, disp_src, &block_trg, data_trg, disp_rank, mirrors_window);
        //     }
    }
    //--------------------------------------------------------------------------
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
void GridBlock::GhostGet_Wait(const Field* field, const lda_t ida, const Wavelet* interp) {
    //--------------------------------------------------------------------------
    const bool do_coarse = (local_parent_.size() + ghost_parent_.size()) > 0;
    if (do_coarse) {
        //................................................
        // copy myself to the coarse block
        {
            // get the neighbor info = source
            const MemSpan      span_src = this->BlockSpan();
            const ConstMemData data_src = this->data(field, ida);

            // get the coarse span = target
            const MemLayout layout_coarse = this->CoarseLayout(interp);
            const MemSpan   span_coarse   = this->CoarseSpan();
            const MemData   data_coarse(&coarse_ptr_, &layout_coarse);

            // copy from me to the coarse block
            const lid_t    shift[3] = {0, 0, 0};
            interp->Copy(1,shift,&span_src,&data_src,&span_coarse,&data_coarse);
            
            // // m_log("reset the coarse block to %d and stride %d ", interp->CoarseNGhostFront(), interp->CoarseStride());
            // const SubBlock coarse_block(interp->CoarseNGhostFront(ghost_len_[0]), interp->CoarseStride(ghost_len_), 0, M_NHALF);
            // // copy myself to the coarse, one point out of 2
            
            // const data_ptr data_src = data(field, ida);
            // const data_ptr data_trg = coarse_ptr_(0, &coarse_block);
            // // interpolate
            // interp->Copy(1, shift, this, data_src, &coarse_block, data_trg);
        }

        //................................................
        // fill the boundary conditions to completely fill the coarse block before the interpolation
        for (auto* const gblock : phys_) {
            // get the direction and the corresponding bctype
            const bctype_t bctype = field->bctype(ida, gblock->iface());
            // // in the face direction, the start and the end are already correct, only the fstart changes
            // SubBlock coarse_block;
            // // interp->CoarseFromFine(gblock, &coarse_block);
            // gblock->ToCoarse(interp, &coarse_block);
            // lid_t fstart[3];
            // // interp->CoarseFromFine(face_start[gblock->iface()], fstart);
            // gblock->ToCoarse(interp, face_start[gblock->iface()], fstart);
            // data_ptr data_trg = coarse_ptr_(0, &coarse_block);

            // get the coarse region
            // const MemLayout layout_coarse = ;
            MemSpan         span_coarse;
            const MemLayout layout_block  = this->BlockLayout();
            const MemLayout layout_coarse = this->CoarseLayout(interp);
            const MemData   data_coarse(&coarse_ptr_, &layout_coarse);
            gblock->GetCoarseSpan(&layout_block, &layout_coarse, &span_coarse);

            // get the start index of the face
            bidx_t fstart[3];
            gblock->GetCoarseLength(&layout_block, &layout_coarse, face_start[gblock->iface()], fstart);

            // apply the BC
            // m_log("apply bc on coarse for block @ %f %f %f", xyz(0), xyz(1), xyz(2));
            if (bctype == M_BC_NEU) {
                NeumanBoundary<M_WAVELET_N - 1> bc;
                bc(gblock->iface(), fstart, hgrid_, 0.0, &span_coarse,& data_coarse);
            } else if (bctype == M_BC_DIR) {
                DirichletBoundary<M_WAVELET_N - 1> bc;
                bc(gblock->iface(), fstart, hgrid_, 0.0,& span_coarse,& data_coarse);
            } else if (bctype == M_BC_EXTRAP) {
                ExtrapBoundary<M_WAVELET_N> bc;
                bc(gblock->iface(), fstart, hgrid_, 0.0, &span_coarse, &data_coarse);
            } else if (bctype == M_BC_ZERO) {
                ZeroBoundary bc;
                bc(gblock->iface(), fstart, hgrid_, 0.0, &span_coarse, &data_coarse);
            } else {
                m_assert(false, "this type of BC is not implemented yet or not valid %d for field <%s>", bctype, field->name().c_str());
            }
        }
        //................................................
        // refine from the coarse tmp to my ghost regions corresponding to coarser neighbors
        {
            // get the coarse span = source
            const MemLayout    layout_coarse = this->CoarseLayout(interp);
            const MemSpan      span_coarse   = this->CoarseExtendedSpan(interp);
            const ConstMemData data_coarse(&coarse_ptr_, &layout_coarse);
            // fixed target quantities
            const MemData data_trg = this->data(field, ida);

            auto coarse_2_parents = [=, &interp](const GhostBlock* gblock) -> void {
                // get the neighbor info = source
                const MemSpan span_trg = gblock->GetSpan();
                // copy from me to the coarse block
                const lid_t shift[3] = {0, 0, 0};
                interp->Interpolate(-1, shift, &span_coarse, &data_coarse, &span_trg, &data_trg);
            };

            for (const auto* gblock : local_parent_) {
                coarse_2_parents(gblock);
            }
            for (const auto* gblock : ghost_parent_) {
                coarse_2_parents(gblock);
            }

            // // take the full coarse block and set the info in my GP
            // const SubBlock block_src(interp->CoarseNGhostFront(ghost_len_[0]),
            //                          interp->CoarseStride(ghost_len_),
            //                          -interp->CoarseNGhostFront(ghost_len_[0]),
            //                          M_NHALF + interp->CoarseNGhostBack(ghost_len_[1]));
            // const data_ptr data_src = coarse_ptr_(0, &block_src);
            // lid_t          shift[3] = {0, 0, 0};

            // for (auto* const gblock : local_parent_) {
            //     // extension from the ghost
            //     // gblock->ExtendGhost(interp->ndetail_citerion_extend_front(), interp->ndetail_citerion_extend_back(), &block_trg);
            //     // interp->Interpolate(-1, shift, &block_src, data_src, &block_trg, data(field, ida));
            //     // no extension
            //     // m_log("@ %e %e %e refinement from %d %d %d to %d %d %d", xyz_[0], xyz_[1], xyz_[2], gblock->start(0), gblock->start(1), gblock->start(2), gblock->end(0), gblock->end(1), gblock->end(2));
            //     interp->Interpolate(-1, shift, &block_src, data_src, gblock, data(field, ida));
            // }
            // for (auto* const gblock : ghost_parent_) {
            //     interp->Interpolate(-1, shift, &block_src, data_src, gblock, data(field, ida));
            // }
        }
    }
    //--------------------------------------------------------------------------
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
void GridBlock::GhostPut_Post(const Field* field, const lda_t ida, const Wavelet* interp, MPI_Win mirrors_window) {
    //--------------------------------------------------------------------------
    const bool do_coarse = (local_parent_.size() + ghost_parent_.size()) > 0;
    if (do_coarse) {
        //................................................
        // apply the physics BC to the best of my knowledge
        // some ghosts are missing but this is okay
        const MemData data_trg = data(field, ida);

        // m_log("apply bc for block @ %f %f %f", xyz(0), xyz(1), xyz(2));
        for (auto* const gblock : phys_) {
            const MemSpan span_bc = gblock->GetSpan();
            bctype_t bctype = field->bctype(ida, gblock->iface());
            // get the correct face_start
            if (bctype == M_BC_NEU) {
                NeumanBoundary<M_WAVELET_N - 1> bc;
                bc(gblock->iface(), face_start[gblock->iface()], hgrid_, 0.0, &span_bc, &data_trg);
            } else if (bctype == M_BC_DIR) {
                DirichletBoundary<M_WAVELET_N - 1> bc;
                bc(gblock->iface(), face_start[gblock->iface()], hgrid_, 0.0,& span_bc, &data_trg);
            } else if (bctype == M_BC_EXTRAP) {
                ExtrapBoundary<M_WAVELET_N> bc;
                bc(gblock->iface(), face_start[gblock->iface()], hgrid_, 0.0,& span_bc,& data_trg);
            } else if (bctype == M_BC_ZERO) {
                ZeroBoundary bc;
                bc(gblock->iface(), face_start[gblock->iface()], hgrid_, 0.0, &span_bc, &data_trg);
            } else {
                m_assert(false, "this type of BC is not implemented yet or not valid %d", bctype);
            }
        }

        //................................................
        // enforce the 0-detail policy in the needed region
        {
            coarse_ptr_.MemSetZero();
            const MemLayout layout = BlockLayout();
            MemData mask(&coarse_ptr_,&layout);
            MemSpan block_span = this->BlockSpan();
            // memset(coarse_ptr_(), 0, CartBlockMemNum(1) * sizeof(real_t));
            // data_ptr mask      = coarse_ptr_(0, this);
            // real_t*  mask_data = mask.Write();

            // lambda to obtain the detail checking pattern
            auto overwrite = [=, &interp](const iface_t ibidule) -> void {
                // get the sign of the ibidule
                real_t sign[3];
                GhostGetSign(ibidule, sign);

                // create the start and end indexes from the current block
                MemSpan smooth_span = block_span;
                // bidx_t smooth_start[3], smooth_end[3];
                for (lda_t id = 0; id < 3; ++id) {
                    if (sign[id] > 0.5) {
                        // my ngh assumed 0 details in my block
                        smooth_span.start[id] = block_span.end[id] -interp->ndetail_citerion_extend_front(); 
                        // smooth_start[ida] = this->end(ida) - interp->ndetail_citerion_extend_front();
                        // the number of my ngh details influencing my values
                        smooth_span.end[id] = block_span.end[id];
                        // smooth_end[ida] = this->end(ida);
                    } else if (sign[id] < (-0.5)) {
                        // my ngh assumed 0 details in my block
                        // smooth_start[ida] = this->start(ida);
                        smooth_span.start[id] = block_span.start[id]; 
                        // the number of my ngh details influencing my values
                        // smooth_end[ida] = this->start(ida) + interp->ndetail_citerion_extend_back();
                        smooth_span.end[id] = block_span.start[id] + interp->ndetail_citerion_extend_back();
                    }
                    // else {
                    //     // even in the directions orthogonal to ibidule, the details must be killed!
                    //     // as my neighbor, which might be fine will kill them as well
                    //     smooth_start[ida] = this->start(ida);
                    //     smooth_end[ida]   = this->end(ida);
                    // }
                }
                // apply it
                // for_loop(&set_mask_to_one, smooth_start, smooth_end);
                // compute the detail, store them in the mask
                // SubBlock block_src(this->gs(), this->stride(), -interp->nghost_front(), M_N + interp->nghost_back());
                // SubBlock block_src(this->gs(), this->stride(), -ghost_len_[0], M_N + ghost_len_[1]);
                // SubBlock block_trg(this->gs(), this->stride(), smooth_start, smooth_end);
                m_assert(ghost_len_[0] >= interp->nghost_front_overwrite(), "the ghost length does not support overwrite: %d vs %d", ghost_len_[0], interp->nghost_front_overwrite());
                m_assert(ghost_len_[1] >= interp->nghost_back_overwrite(), "the ghost length does not support overwrite: %d vs %d", ghost_len_[1], interp->nghost_back_overwrite());
                const MemSpan me       = ExtendedSpan(ghost_len_);
                const MemData data_trg = this->data(field, ida);
                interp->OverwriteDetails(&me, &smooth_span, &data_trg);
                // interp->OverwriteDetails(&block_src, &block_trg, this->data(field, ida));
            };
            // let's god
            for (auto* gblock : local_parent_) {
                overwrite(gblock->ibidule());
            }
            for (auto* gblock : ghost_parent_) {
                overwrite(gblock->ibidule());
            }
        }
        //................................................
        // coarsen to put to my finer neighbors
        {
            //................................................
            // reset the tmp to use for the put operations
            coarse_ptr_.MemSetZero();
            // memset(coarse_ptr_(), 0, interp->CoarseSize(ghost_len_) * sizeof(real_t));

            // // I am now complete (except children GP), get my coarse representation
            const MemLayout coarse_layout = this->CoarseLayout(interp);
            const MemSpan   span_coarse   = this->CoarseSpan();
            const MemData   data_coarse(&coarse_ptr_, &coarse_layout);
            // const SubBlock coarse_block(interp->CoarseNGhostFront(ghost_len_[0]), interp->CoarseStride(ghost_len_), 0, M_NHALF);
            // data_ptr       data_coarse = coarse_ptr_(0, &coarse_block);

            // the source block is the ghost extended block
            const lid_t shift[3] = {0, 0, 0};
            // const SubBlock me_extended(M_GS, M_STRIDE, - ghost_len_[0], M_N + ghost_len_[1]);
            const MemSpan ext_span = this->ExtendedSpan(ghost_len_);
            // interpolate, the level is 1 coarser and the shift is unchanged
            const ConstMemData data_src = this->ConstData(field, ida);
            interp->Interpolate(1, shift, &ext_span, &data_src, &span_coarse, &data_coarse);
            // interp->Interpolate(1, shift, &me_extended, data(field, ida), &coarse_block, data_coarse);
            m_assert(ghost_len_[0] >= interp->nghost_front_coarsen() && ghost_len_[0] >= interp->nghost_back_coarsen(), "the ghost sizes %d %d must be >= %d %d", ghost_len_[0], ghost_len_[1], interp->nghost_front_coarsen(), interp->nghost_back_coarsen());

            //................................................
            // start the put commands to set the value to my neighbors
            // loop on the ghost list
            for (auto* const gblock : ghost_parent_reverse_) {
                const MPI_Aint  disp_trg   = gblock->data_src();
                const rank_t    trg_rank   = gblock->rank();
                const MemLayout layout_trg = gblock->SourceLayout();
                const MemSpan   span_trg   = gblock->GetSpan();
                // interpolate, the parent's mirror have been created to act on the tmp
                m_assert(gblock->dlvl() == 0, "we must have a level 0, here %d", gblock->dlvl());
                const ConstMemData data_src(&coarse_ptr_, &coarse_layout);
                PutRma(gblock->dlvl(), gblock->shift(), &coarse_layout, &span_coarse, &data_src, &layout_trg, &span_trg, disp_trg, trg_rank, mirrors_window);
                // interp->PutRma(gblock->dlvl(), gblock->shift(), &coarse_block, data_coarse, gblock, disp_trg, trg_rank, mirrors_window);
            }
            for (auto* const gblock : local_parent_reverse_) {
                const CartBlock* ngh_block = gblock->data_src();
                const MemSpan    span_trg  = gblock->GetSpan();
                const MemData    data_trg  = ngh_block->data(field, ida);
                // interpolate, the level is 1 coarser and the shift is unchanged
                m_assert(gblock->dlvl() == 0, "we must have a level 0, here %d", gblock->dlvl());
                const ConstMemData data_src(&coarse_ptr_, &coarse_layout);
                interp->Copy(gblock->dlvl(), gblock->shift(), &span_coarse, &data_src, &span_trg, &data_trg);
                // interp->Copy(gblock->dlvl(), gblock->shift(), &coarse_block, data_coarse, gblock, data_trg);
            }
        }
    }
    //--------------------------------------------------------------------------
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
void GridBlock::GhostPut_Wait(const Field* field, const lda_t ida, const Wavelet* interp) {
    //-------------------------------------------------------------------------
    const MemData data_trg = data(field, ida);
    // data_ptr data_trg = data(field, ida);
    // m_log("apply bc for block @ %f %f %f", xyz(0), xyz(1), xyz(2));
    for (auto gblock : phys_) {
        const MemSpan span_bc = gblock->GetSpan();
        bctype_t bctype = field->bctype(ida, gblock->iface());
        // get the correct face_start
        if (bctype == M_BC_NEU) {
            NeumanBoundary<M_WAVELET_N - 1> bc;
            bc(gblock->iface(), face_start[gblock->iface()], hgrid_, 0.0, &span_bc,& data_trg);
        } else if (bctype == M_BC_DIR) {
            DirichletBoundary<M_WAVELET_N - 1> bc;
            bc(gblock->iface(), face_start[gblock->iface()], hgrid_, 0.0, &span_bc, &data_trg);
        } else if (bctype == M_BC_EXTRAP) {
            ExtrapBoundary<M_WAVELET_N> bc;
            bc(gblock->iface(), face_start[gblock->iface()], hgrid_, 0.0, &span_bc, &data_trg);
        } else if (bctype == M_BC_ZERO) {
            ZeroBoundary bc;
            bc(gblock->iface(), face_start[gblock->iface()], hgrid_, 0.0, &span_bc, &data_trg);
        } else {
            m_assert(false, "this type of BC is not implemented yet or not valid %d", bctype);
        }
        //--------------------------------------------------------------------------
    }
    //--------------------------------------------------------------------------
}


// /**
//  * @brief Downsample the block and re-evaluate the boundary conditions on the coarse version
//  *
//  * @param field
//  * @param ida
//  * @param interp
//  */
// void GridBlock::Coarse_DownSampleWithBoundary(const Field* field, const lda_t ida, const Wavelet* interp, SubBlock* coarse_block) {
//     //--------------------------------------------------------------------------
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
//     //--------------------------------------------------------------------------
// }
