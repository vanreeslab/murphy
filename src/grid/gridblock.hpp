#ifndef SRC_GRID_GRIDBLOCK_HPP_
#define SRC_GRID_GRIDBLOCK_HPP_

#include <p8est.h>

#include <list>

#include "core/macros.hpp"
#include "core/pointers.hpp"
#include "core/types.hpp"
#include "grid/cartblock.hpp"
#include "grid/forestgrid.hpp"
#include "grid/ghostblock.hpp"
#include "grid/physblock.hpp"
#include "tools/patch.hpp"
#include "tools/prof.hpp"
#include "wavelet/wavelet.hpp"

#define M_NNEIGHBORS 26

typedef enum StatusAdapt {
    M_ADAPT_NONE,
    M_ADAPT_COARSER,
    M_ADAPT_FINER,
    M_ADAPT_SAME,
    M_ADAPT_NEW_COARSE,
    M_ADAPT_NEW_FINE
} StatusAdapt;

/**
 * @brief a @ref CartBlock with ghosting and wavelet capabilities
 * 
 */
class GridBlock : public CartBlock {
   protected:
    StatusAdapt status_lvl_ = M_ADAPT_NONE;   //!< indicat the status of the level: +1 should refine, -1 should coarsen
    short_t*    status_ngh_;                  //!< indicate if my sibling neighbors are going to coarsen, stored as local_sibling and then ghost_sibling
    short_t     n_dependency_active_ = 0;     //!< list of dependency = how to create my information after refinement/coarsening
    GridBlock*  dependency_[P8EST_CHILDREN];  //!< the pointer to the dependency block

    mem_ptr coarse_ptr_;  //!< a memory reserved for coarser version of myself, includes ghost points

    // list of ghosting
    std::list<GhostBlock<GridBlock*>*> local_sibling_;         //<! local neighbors at my resolution
    std::list<GhostBlock<GridBlock*>*> local_parent_;          //!< local neighbors coarser (neighbor to me)
    std::list<GhostBlock<GridBlock*>*> local_children_;        //!< local neighbors finer
    std::list<GhostBlock<GridBlock*>*> local_parent_reverse_;  //!< local neighbors coarser (me to neighbors)
    std::list<GhostBlock<MPI_Aint>*>   ghost_sibling_;         //<! ghost neighbors at my resolution
    std::list<GhostBlock<MPI_Aint>*>   ghost_parent_;          //!< ghost neighbors coarser (neighbor to me)
    std::list<GhostBlock<MPI_Aint>*>   ghost_children_;        //!< ghost neighbors coarser (neighbor to me)
    std::list<GhostBlock<MPI_Aint>*>   ghost_parent_reverse_;  //!<  ghost neighbors coarser (me to neighbors)
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
    void        ResetStatus() { status_lvl_ = M_ADAPT_NONE; };

    void UpdateStatusFromCriterion(/* params */ m_ptr<const Wavelet> interp, const real_t rtol, const real_t ctol, m_ptr<const Field> field_citerion,
                                   /* prof */ m_ptr<Prof> profiler);
    void UpdateStatusFromPatches(/* params */ m_ptr<const Wavelet> interp, m_ptr<std::list<Patch> > patch_list,
                                 /* prof */ m_ptr<Prof> profiler);

    void ComputeDetails(m_ptr<const Wavelet> interp, m_ptr<const Field> criterion, m_ptr<const Field> details);
    void UpdateSmoothingMask(const m_ptr<const Wavelet>& interp);

    void FWTAndGetStatus(m_ptr<const Wavelet> interp, const real_t rtol, const real_t ctol, m_ptr<const Field> field_citerion, m_ptr<Prof> profiler);
    void SetNewByCoarsening(m_ptr<const qid_t> qid, const m_ptr<short_t> coarsen_vec) const;
    void GetNewByCoarseningFromNeighbors(const m_ptr<const short_t> status_vec, MPI_Win status_window);
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
    void       SolveDependency(m_ptr<const Wavelet> interp, std::map<std::string, m_ptr<Field> >::const_iterator field_start, std::map<std::string, m_ptr<Field> >::const_iterator field_end, m_ptr<Prof> profiler);
    void       SmoothResolutionJump(m_ptr<const Wavelet> interp, std::map<std::string, m_ptr<Field> >::const_iterator field_start, std::map<std::string, m_ptr<Field> >::const_iterator field_end, m_ptr<Prof> profiler);
    void       ClearResolutionJump(m_ptr<const Wavelet> interp, std::map<std::string, m_ptr<Field> >::const_iterator field_start, std::map<std::string, m_ptr<Field> >::const_iterator field_end, m_ptr<Prof> profiler);
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
    std::list<GhostBlock<GridBlock*>*>* local_sibling() { return &local_sibling_; };
    std::list<GhostBlock<GridBlock*>*>* local_parent() { return &local_parent_; };
    std::list<GhostBlock<GridBlock*>*>* local_children() { return &local_children_; };
    std::list<GhostBlock<GridBlock*>*>* local_parent_reverse() { return &local_parent_reverse_; };
    std::list<GhostBlock<MPI_Aint>*>*   ghost_sibling() { return &ghost_sibling_; };
    std::list<GhostBlock<MPI_Aint>*>*   ghost_parent() { return &ghost_parent_; };
    std::list<GhostBlock<MPI_Aint>*>*   ghost_children() { return &ghost_children_; };
    std::list<GhostBlock<MPI_Aint>*>*   ghost_parent_reverse() { return &ghost_parent_reverse_; };
    std::list<PhysBlock*>*              phys() { return &phys_; };
    /**@} */

    void GhostInitLists(m_ptr<const qid_t> qid, m_ptr<const ForestGrid> grid, m_ptr<const Wavelet> interp, MPI_Win local2disp_window);
    void GhostGet_Cmpt(m_ptr<const Field> field, const lda_t ida, m_ptr<const Wavelet> interp);
    void GhostGet_Post(m_ptr<const Field> field, const lda_t ida, m_ptr<const Wavelet> interp, MPI_Win mirrors_window);
    void GhostGet_Wait(m_ptr<const Field> field, const lda_t ida, m_ptr<const Wavelet> interp);
    void GhostPut_Post(m_ptr<const Field> field, const lda_t ida, m_ptr<const Wavelet> interp, MPI_Win mirrors_window);
    void GhostPut_Wait(m_ptr<const Field> field, const lda_t ida, m_ptr<const Wavelet> interp);
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
