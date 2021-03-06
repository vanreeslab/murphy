#ifndef SRC_GRID_GRIDBLOCK_HPP_
#define SRC_GRID_GRIDBLOCK_HPP_

#include <p8est.h>

#include <list>

#include "core/macros.hpp"
#include "core/types.hpp"
#include "grid/cartblock.hpp"
#include "grid/neighborblock.hpp"
#include "grid/physblock.hpp"
#include "tools/patch.hpp"
#include "tools/toolsp4est.hpp"
#include "wavelet/wavelet.hpp"

#define M_NNEIGHBORS 26

typedef enum StatusAdapt {
    M_ADAPT_NEW_COARSE,
    M_ADAPT_NEW_FINE,
    M_ADAPT_SAME,
    M_ADAPT_COARSER,
    M_ADAPT_FINER,
} StatusAdapt;

using GBLocal      = NeighborBlock<CartBlock*>;
using GBMirror     = NeighborBlock<MPI_Aint>;
using GBPhysic     = PhysBlock;
using ListGBLocal  = std::list<GBLocal*>;
using ListGBMirror = std::list<GBMirror*>;
using ListGBPhysic = std::list<GBPhysic*>;

typedef enum StatusNghIndex {
    M_LOC_PARENT   = 0,
    M_GLO_PARENT   = 1,
    M_LOC_SIBLING  = 2,
    M_GLO_SIBLING  = 3,
    M_LOC_CHILDREN = 4,
    M_GLO_CHILDREN = 5,

} StatusNghIndex;


/**
 * @brief a @ref CartBlock with ghosting and wavelet capabilities
 * 
 */
class GridBlock : public CartBlock {
   private:
    // active ghost lengths
    bidx_t ghost_len_[2] = {0, 0};  //!< contains the current ghost length

    // status tracking
    bool        status_refined_ = false;         //!< track if the block has been refined, which prevents a new coarsening
    StatusAdapt status_lvl_     = M_ADAPT_SAME;  //!< indicate the status of the block
    short_t*    status_ngh_[6]  = {nullptr};     //!< status of my neighbors, see StatusNghIndex

    // dependency tracking
    short_t    n_dependency_active_        = 0;          //!< list of dependency = how to create my information after refinement/coarsening
    GridBlock* dependency_[P8EST_CHILDREN] = {nullptr};  //!< the pointer to the dependency block

    MemPtr coarse_ptr_;  //!< a memory reserved for coarser version of myself, includes ghost points

    // list of ghosting
    ListGBLocal  local_sibling_;         //<! local neighbors at my resolution
    ListGBLocal  local_parent_;          //!< local neighbors coarser (neighbor to me)
    ListGBLocal  local_children_;        //!< local neighbors finer
    ListGBLocal  local_parent_reverse_;  //!< local neighbors coarser (me gto neighbors)
    ListGBMirror ghost_sibling_;         //<! ghost neighbors at my resolution
    ListGBMirror ghost_parent_;          //!< ghost neighbors coarser (neighbor to me)
    ListGBMirror ghost_children_;        //!< ghost neighbors coarser (neighbor to me)
    ListGBMirror ghost_parent_reverse_;  //!<  ghost neighbors coarser (me to neighbors)
    ListGBPhysic phys_;                  //!<  physical boundary condition

   public:
    explicit GridBlock(const real_t length, const real_t xyz[3], const sid_t level);
    virtual ~GridBlock();

    M_INLINE MemLayout CoarseLayout(const Wavelet* interp) const {
        return MemLayout(M_LAYOUT_BLOCK, interp->CoarseNGhostFront(ghost_len_[0]), M_NHALF, interp->CoarseNGhostBack(ghost_len_[1]));
    }
    M_INLINE MemSpan CoarseSpan() const {
        return MemSpan(0, M_NHALF);
    }

    M_INLINE MemSpan CoarseExtendedSpan(const Wavelet* interp) const {
        return MemSpan(-interp->CoarseNGhostFront(ghost_len_[0]), M_NHALF + interp->CoarseNGhostBack(ghost_len_[1]));
    }

    // M_INLINE MemSpan ExtendedSpanCurrent() const {
    //     return MemSpan(-ghost_len_[0], M_N + ghost_len_[1]);
    // }
    M_INLINE MemSpan ExtendedSpan(const bidx_t ghost_len[2]) const {
        return MemSpan(-ghost_len[0], M_N + ghost_len[1]);
    }

    /**
     * @brief return the offset for the partitioning buffer
     * 
     * @warning the offset is measured in number of real_t
     */
    virtual size_t PartitionDataOffset() const override { return 2; }

    /**
     * @brief pack (store) the needed information in the partitioner buffer
     * 
     * The data will be sent over to another rank during partitioning and recoverd by PartitionDataUnPack()
     * 
     * @warning the information MUST be cast to a real_t
     */
    virtual void PartitionDataPack(real_t* buff) const override {
        m_assert(sizeof(status_lvl_) < sizeof(real_t), "the size of the status must fit in the real type");
        m_assert(sizeof(status_refined_) < sizeof(real_t), "the size of the status must fit in the real type");
        buff[0] = static_cast<real_t>(status_lvl_);
        buff[1] = static_cast<real_t>(status_refined_);
    }

    /**
     * @brief unpack the needed information from the partitioner buffer
     * 
     * The data has been stored by PartitionDataPack() and send over to another rank
     * 
     * @warning the information MUST be uncast from a real_t
     */
    virtual void PartitionDataUnPack(const real_t* buff) override {
        m_assert(sizeof(status_lvl_) < sizeof(real_t), "the size of the status must fit in the real type");
        m_assert(sizeof(status_refined_) < sizeof(real_t), "the size of the status must fit in the real type");
        status_lvl_     = static_cast<StatusAdapt>(buff[0]);
        status_refined_ = static_cast<bool>(buff[1]);
    }

    /**
     * @name Status level management
     *
     * @{ */

    void status_level(const StatusAdapt status) { status_lvl_ = status; };
    void status_refined(const bool status) { status_refined_ = status; }

    [[nodiscard]] StatusAdapt status_level() const { return status_lvl_; };
    [[nodiscard]] bool        status_refined() const { return status_refined_; }

    /** @brief reset the status, ready to go for adaptation */
    void StatusReset() {
        m_assert(M_ADAPT_NEW_FINE < M_ADAPT_SAME && M_ADAPT_NEW_COARSE < M_ADAPT_SAME, "please keep M_ADAPT_NEW_FINE/COARSE/M_ADAPT_NONE < M_ADAPT_SAME");
        status_lvl_     = M_ADAPT_SAME;
        status_refined_ = false;
    };

    /** @brief reset the status after one pass of adaptation, register if the block has been refined already once */
    void StatusCleanup() {
        m_assert(M_ADAPT_NEW_FINE < M_ADAPT_SAME && M_ADAPT_NEW_COARSE < M_ADAPT_SAME, "please keep M_ADAPT_NEW_FINE/COARSE/M_ADAPT_NONE < M_ADAPT_SAME");
        status_refined_ = status_refined_ || (status_lvl_ == M_ADAPT_NEW_FINE);
        status_lvl_     = M_ADAPT_SAME;
    };


    void UpdateStatusFromCriterion(const Wavelet* interp, const real_t rtol, const real_t ctol, const Field* field_citerion, const lda_t ida);
    void UpdateStatusFromPatches(const Wavelet* interp, std::list<Patch>* patch_list);
    
    void UpdateStatusFromLevel(const level_t min_level, const level_t max_level);
    void UpdateStatusForwardRefinement();
    void UpdateStatusFromGlobalPolicy();

    void MaxMinDetails(const Wavelet* interp, const Field* criterion, real_t maxmin[2],
                       bidx_t* max_blocks, const real_t max_cat, const real_t min_cat, const short_t n_cat);
    void StoreDetails(const Wavelet* interp, const Field* criterion, const Field* details);
    
    void SyncStatusInit();
    void SyncStatusFill(const qid_t* qid, short_t* const coarsen_vec);
    void SyncStatusUpdate(const short_t* const status_vec, MPI_Win status_window);
    void SyncStatusFinalize();
    /** @} */

    /**
     * @name dependency management
     * 
     * @{
     */
    sid_t      n_dependency_active() { return n_dependency_active_; }
    GridBlock* PopDependency(const sid_t child_id);
    void       PushDependency(const sid_t child_id, GridBlock* dependent_block);
    void       SolveDependency(const Wavelet* interp, std::map<std::string, Field*>::const_iterator field_start, std::map<std::string, Field*>::const_iterator field_end);  //, Prof*  profiler);
    void       SmoothResolutionJump(const Wavelet* interp, std::map<std::string, Field*>::const_iterator field_start,
                                    std::map<std::string, Field*>::const_iterator field_end);
    /** @} */

    /**
     * @name return the ghost lists
     * @{
     */
    ListGBLocal*  local_sibling() { return &local_sibling_; };
    ListGBLocal*  local_parent() { return &local_parent_; };
    ListGBLocal*  local_children() { return &local_children_; };
    ListGBLocal*  local_parent_reverse() { return &local_parent_reverse_; };
    ListGBMirror* ghost_sibling() { return &ghost_sibling_; };
    ListGBMirror* ghost_parent() { return &ghost_parent_; };
    ListGBMirror* ghost_children() { return &ghost_children_; };
    ListGBMirror* ghost_parent_reverse() { return &ghost_parent_reverse_; };
    ListGBPhysic* phys() { return &phys_; };
    /**@} */

    void GhostUpdateSize(const bidx_t ghost_len[2]);

    void GhostInitLists(const qid_t* qid, const p4est_Essentials* ess_info, const Wavelet* interp, MPI_Win local2disp_window);
    void GhostGet_Cmpt(const Field* field, const lda_t ida, const Wavelet* interp);
    void GhostGet_Post(const Field* field, const lda_t ida, const Wavelet* interp, MPI_Win mirrors_window);
    void GhostGet_Wait(const Field* field, const lda_t ida, const Wavelet* interp);
    void GhostPut_Post(const Field* field, const lda_t ida, const Wavelet* interp, MPI_Win mirrors_window);
    void GhostPut_Wait(const Field* field, const lda_t ida, const Wavelet* interp);
    void GhostFreeLists();
};

#endif  // SRC_GRIDBLOCK_HPP_
