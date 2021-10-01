#ifndef SRC_GRID_GRIDBLOCK_HPP_
#define SRC_GRID_GRIDBLOCK_HPP_

#include <p8est.h>

#include <list>

#include "core/macros.hpp"
#include "core/types.hpp"
#include "grid/cartblock.hpp"
#include "grid/forestgrid.hpp"
#include "grid/neighborblock.hpp"
#include "grid/physblock.hpp"
#include "tools/patch.hpp"
#include "tools/prof.hpp"
#include "wavelet/wavelet.hpp"

#define M_NNEIGHBORS 26

typedef enum StatusAdapt {
    M_ADAPT_NEW_COARSE,
    M_ADAPT_NEW_FINE,
    M_ADAPT_SAME,
    M_ADAPT_COARSER,
    M_ADAPT_FINER,
} StatusAdapt;

// using GBLocal      = NeighborBlock<GridBlock*>;
using GBLocal      = NeighborBlock<CartBlock*>;
using GBMirror     = NeighborBlock<MPI_Aint>;
using GBPhysic     = PhysBlock;
using ListGBLocal  = std::list<GBLocal*>;
using ListGBMirror = std::list<GBMirror*>;
using ListGBPhysic = std::list<GBPhysic*>;

typedef enum StatusNghIndex {
    M_LOC_PARENT = 0,
    M_GLO_PARENT = 1,
    M_LOC_SIBLING = 2,
    M_GLO_SIBLING = 3,
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
    ~GridBlock();

    __attribute__((always_inline)) inline MemLayout CoarseLayout(const Wavelet* interp) const {
        return MemLayout(M_LAYOUT_BLOCK, interp->CoarseNGhostFront(ghost_len_[0]), M_NHALF, interp->CoarseNGhostBack(ghost_len_[1]));
    }
    __attribute__((always_inline)) inline MemSpan CoarseSpan() const {
        return MemSpan(0, M_NHALF);
    }

    __attribute__((always_inline)) inline MemSpan CoarseExtendedSpan(const Wavelet* interp) const {
        return MemSpan(-interp->CoarseNGhostFront(ghost_len_[0]), M_NHALF + interp->CoarseNGhostBack(ghost_len_[1]));
    }

    // __attribute__((always_inline)) inline MemSpan ExtendedSpanCurrent() const {
    //     return MemSpan(-ghost_len_[0], M_N + ghost_len_[1]);
    // }
    __attribute__((always_inline)) inline MemSpan ExtendedSpan(const bidx_t ghost_len[2]) const {
        return MemSpan(-ghost_len[0], M_N + ghost_len[1]);
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

    // void StatusCleanup() { status_lvl_ = (status_lvl_ == M_ADAPT_NONE)? M_ADAPT_SAME : status_lvl_; };

    /** @brief set the status to M_ADAPT_NONE unless it is M_ADAPT_NEW_FINE/COARSE */
    // void StatusRememberPast() {
    //     // will preserve the information "newly" refined/coarsed
        
    //     status_lvl_ = m_min(M_ADAPT_NONE, status_lvl_);
    // };
    // /** @brief set the status to M_ADAPT_SAME if the status is from the past (i.e. M_ADAPT_NEW_FINE/COARSE) or if we still have M_ADAPT_NONE  */
    // void StatusForgetPast() {
    //     m_assert(M_ADAPT_NEW_FINE < M_ADAPT_SAME && M_ADAPT_NEW_COARSE < M_ADAPT_SAME, "please keep M_ADAPT_NEW_FINE/COARSE < M_ADAPT_SAME");
    //     m_assert(M_ADAPT_SAME < M_ADAPT_COARSER && M_ADAPT_SAME < M_ADAPT_FINER, "please keep M_ADAPT_SAME < M_ADAPT_FINER/COARSER");
    //     m_assert(M_ADAPT_NONE < M_ADAPT_SAME && M_ADAPT_NONE < M_ADAPT_FINER && M_ADAPT_NONE < M_ADAPT_COARSER, "please keep M_ADAPT_NONE < M_ADAPT_SAME < M_ADAPT_FINER/COARSER");
    //     status_lvl_ = m_max(M_ADAPT_SAME, status_lvl_);
    // };

    void UpdateStatusFromCriterion(const Wavelet* interp, const real_t rtol, const real_t ctol, const Field* field_citerion, const lda_t ida);
    void UpdateStatusFromPatches(const Wavelet* interp, std::list<Patch>* patch_list);
    // void UpdateSmoothingMask(const Wavelet* interp);
    
    void UpdateStatusFromLevel(const level_t min_level, const level_t max_level);
    void UpdateStatusForwardRefinement();
    void UpdateStatusFromGlobalPolicy();

    void MaxMinDetails(const Wavelet* interp, const Field* criterion, real_t maxmin[2],
                       bidx_t* max_blocks, const real_t max_cat, const real_t min_cat, const short_t n_cat);
    void StoreDetails(const Wavelet* interp, const Field* criterion, const Field* details);
    // void UpdateSmoothingMask(const Wavelet* const interp);

    // void FWTAndGetStatus(const Wavelet*  interp, const real_t rtol, const real_t ctol, const Field*  field_citerion, Prof*  profiler);
    void SyncStatusInit();
    void SyncStatusFill(const qid_t* qid, short_t* const coarsen_vec);
    void SyncStatusUpdate(const short_t* const status_vec, MPI_Win status_window);
    void SyncStatusFinalize();
    // void UpdateDetails();
    /** @} */

    /**
     * @name dependency management
     * 
     * @{
     */
    sid_t      n_dependency_active() { return n_dependency_active_; }
    GridBlock* PopDependency(const sid_t child_id);
    void       PushDependency(const sid_t child_id, GridBlock* dependent_block);
    void       SolveDependency(const Wavelet*  interp, std::map<std::string, Field*  >::const_iterator field_start, std::map<std::string, Field*  >::const_iterator field_end);//, Prof*  profiler);
    void       SmoothResolutionJump(const Wavelet* interp, std::map<std::string, Field*>::const_iterator field_start, 
                                                            std::map<std::string, Field*>::const_iterator field_end);
    // void       ClearResolutionJump(const Wavelet*  interp, std::map<std::string, Field*  >::const_iterator field_start, std::map<std::string, Field*  >::const_iterator field_end, Prof*  profiler);
    /** @} */

    // /**
    //  * @name coarse pointer management
    //  * @{
    //  */
    // Mem coarse_ptr() const { return coarse_ptr_; }
    // void    AllocateCoarsePtr(const size_t memsize);
    // /**@} */

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

    void GhostInitLists(const qid_t* qid, const ForestGrid* grid, const Wavelet* interp, MPI_Win local2disp_window);
    void GhostGet_Cmpt(const Field* field, const lda_t ida, const Wavelet* interp);
    void GhostGet_Post(const Field* field, const lda_t ida, const Wavelet* interp, MPI_Win mirrors_window);
    void GhostGet_Wait(const Field* field, const lda_t ida, const Wavelet* interp);
    void GhostPut_Post(const Field* field, const lda_t ida, const Wavelet* interp, MPI_Win mirrors_window);
    void GhostPut_Wait(const Field* field, const lda_t ida, const Wavelet* interp);
    void GhostFreeLists();

    // void Coarse_DownSampleWithBoundary(const Field* field, const lda_t ida, const Wavelet* interp, SubBlock* coarse_block);
};

static inline GridBlock* p4est_GetGridBlock(const qdrt_t* quad) {
    // the user_data is an array with the addresses of the block
    GridBlock** p4est_usr_data = static_cast<GridBlock**>(quad->p.user_data);
    GridBlock*  block          = p4est_usr_data[0];
    m_assert(block != nullptr, "the block address cannot be null, we have an issue here");
    return block;
}

static inline void p4est_SetGridBlock(qdrt_t* quad, GridBlock* block) {
    // reinterpret_cast<GridBlock**>(quad->p.user_data)[0] = block;
    m_assert(block != nullptr, "the block address cannot be null, we have an issue here");
    GridBlock** p4est_usr_data = static_cast<GridBlock**>(quad->p.user_data);
    p4est_usr_data[0]          = block;
}

#endif  // SRC_GRIDBLOCK_HPP_
