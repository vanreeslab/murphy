#ifndef SRC_GRID_GRIDBLOCK_HPP_
#define SRC_GRID_GRIDBLOCK_HPP_

#include <p8est.h>

#include <list>

#include "core/macros.hpp"
#include "core/pointers.hpp"
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
    M_ADAPT_NONE,
    M_ADAPT_SAME,
    M_ADAPT_COARSER,
    M_ADAPT_FINER
} StatusAdapt;

/**
 * @brief a @ref CartBlock with ghosting and wavelet capabilities
 * 
 */
class GridBlock : public CartBlock {
   private:
    bidx_t      ghost_len_[2];                               //!< contains the current ghost length
    StatusAdapt status_lvl_                 = M_ADAPT_NONE;  //!< indicate the status of the block
    short_t*    status_ngh_                 = nullptr;       //!< indicate if my sibling neighbors are going to coarsen, stored as local_sibling and then ghost_sibling
    short_t     n_dependency_active_        = 0;             //!< list of dependency = how to create my information after refinement/coarsening
    GridBlock*  dependency_[P8EST_CHILDREN] = {nullptr};     //!< the pointer to the dependency block

    mem_ptr coarse_ptr_;  //!< a memory reserved for coarser version of myself, includes ghost points

    // list of ghosting
    std::list<NeighborBlock<GridBlock*>*> local_sibling_;         //<! local neighbors at my resolution
    std::list<NeighborBlock<GridBlock*>*> local_parent_;          //!< local neighbors coarser (neighbor to me)
    std::list<NeighborBlock<GridBlock*>*> local_children_;        //!< local neighbors finer
    std::list<NeighborBlock<GridBlock*>*> local_parent_reverse_;  //!< local neighbors coarser (me to neighbors)
    std::list<NeighborBlock<MPI_Aint>*>   ghost_sibling_;         //<! ghost neighbors at my resolution
    std::list<NeighborBlock<MPI_Aint>*>   ghost_parent_;          //!< ghost neighbors coarser (neighbor to me)
    std::list<NeighborBlock<MPI_Aint>*>   ghost_children_;        //!< ghost neighbors coarser (neighbor to me)
    std::list<NeighborBlock<MPI_Aint>*>   ghost_parent_reverse_;  //!<  ghost neighbors coarser (me to neighbors)
    std::list<PhysBlock*>              phys_;                  //!<  physical boundary condition

   public:
    explicit GridBlock(const real_t length, const real_t xyz[3], const sid_t level);
    ~GridBlock();

    /**
     * @name Status level management
     *
     * @{ */
    StatusAdapt status_level() const { return status_lvl_; };
    void        status_level(const StatusAdapt status) { status_lvl_ = status; };

    /** @brief set the status to M_ADAPT_NONE*/
    void StatusReset() { status_lvl_ = M_ADAPT_NONE; };

    /** @brief set the status to M_ADAPT_NONE unless it is M_ADAPT_NEW_FINE/COARSE */
    void StatusRememberPast() {
        // will preserve the information "newly" refined/coarsed
        status_lvl_ = m_min(M_ADAPT_NONE, status_lvl_);
    };
    /** @brief set the status to M_ADAPT_SAME if the status is from the past (i.e. M_ADAPT_NEW_FINE/COARSE) or if we still have M_ADAPT_NONE  */
    void StatusForgetPast() {
        status_lvl_ = m_max(M_ADAPT_SAME, status_lvl_);
    };

    void UpdateStatusFromCriterion(/* params */ const lda_t ida, const Wavelet* interp, const real_t rtol, const real_t ctol, const Field* field_citerion,
                                   /* prof */ Prof* profiler);
    void UpdateStatusFromPatches(/* params */ const Wavelet* interp, std::list<Patch>* patch_list,
                                 /* prof */ Prof* profiler);
    void UpdateStatusFromPolicy();

    void StoreDetails(const Wavelet* interp, const Field* criterion, const Field* details);
    void UpdateSmoothingMask(const Wavelet* const interp);

    // void FWTAndGetStatus(const Wavelet*  interp, const real_t rtol, const real_t ctol, const Field*  field_citerion, Prof*  profiler);
    void SetNewByCoarsening(const qid_t*  qid, short_t* const  coarsen_vec) const;
    void GetNewByCoarseningFromNeighbors(const short_t* const  status_vec, MPI_Win status_window);
    void UpdateDetails();
    /** @} */

    /**
     * @name dependency management
     * 
     * @{
     */
    sid_t      n_dependency_active() { return n_dependency_active_; }
    GridBlock* PopDependency(const sid_t child_id);
    void       PushDependency(const sid_t child_id, GridBlock* dependent_block);
    void       SolveDependency(const Wavelet*  interp, std::map<std::string, Field*  >::const_iterator field_start, std::map<std::string, Field*  >::const_iterator field_end, Prof*  profiler);
    void       SmoothResolutionJump(const Wavelet*  interp, std::map<std::string, Field*  >::const_iterator field_start, std::map<std::string, Field*  >::const_iterator field_end, Prof*  profiler);
    // void       ClearResolutionJump(const Wavelet*  interp, std::map<std::string, Field*  >::const_iterator field_start, std::map<std::string, Field*  >::const_iterator field_end, Prof*  profiler);
    /** @} */

    /**
     * @name coarse pointer management
     * @{
     */
    mem_ptr coarse_ptr() const { return coarse_ptr_; }
    void    AllocateCoarsePtr(const size_t memsize);
    /**@} */

    /**
     * @name return the ghost lists
     * @{
     */
    std::list<NeighborBlock<GridBlock*>*>* local_sibling() { return &local_sibling_; };
    std::list<NeighborBlock<GridBlock*>*>* local_parent() { return &local_parent_; };
    std::list<NeighborBlock<GridBlock*>*>* local_children() { return &local_children_; };
    std::list<NeighborBlock<GridBlock*>*>* local_parent_reverse() { return &local_parent_reverse_; };
    std::list<NeighborBlock<MPI_Aint>*>*   ghost_sibling() { return &ghost_sibling_; };
    std::list<NeighborBlock<MPI_Aint>*>*   ghost_parent() { return &ghost_parent_; };
    std::list<NeighborBlock<MPI_Aint>*>*   ghost_children() { return &ghost_children_; };
    std::list<NeighborBlock<MPI_Aint>*>*   ghost_parent_reverse() { return &ghost_parent_reverse_; };
    std::list<PhysBlock*>*                 phys() { return &phys_; };
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
    return *(reinterpret_cast<GridBlock**>(quad->p.user_data));
}

static inline void p4est_SetGridBlock(qdrt_t* quad, GridBlock* block) {
    *(reinterpret_cast<GridBlock**>(quad->p.user_data)) = block;
}

#endif  // SRC_GRIDBLOCK_HPP_
