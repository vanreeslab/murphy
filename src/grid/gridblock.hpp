#ifndef SRC_GRID_GRIDBLOCK_HPP_
#define SRC_GRID_GRIDBLOCK_HPP_

#include <mpi.h>

#include <array>
#include <limits>
#include <list>
#include <string>
#include <unordered_map>

#include "core/macros.hpp"
#include "core/memlayout.hpp"
#include "core/pointers.hpp"
#include "core/types.hpp"
#include "field.hpp"
#include "forestgrid.hpp"
#include "ghostblock.hpp"
#include "p8est.h"
#include "physblock.hpp"
#include "prof.hpp"
#include "wavelet/wavelet.hpp"

/**
 * @brief implements a @ref Block that is used as a leaf for the tree
 * 
 */
class GridBlock : public MemLayout {
   protected:
    sid_t   status_lvl_ = 0;                //!< indicate if the block has to change: +1 -> must be refined, -1 must be coarsened, 0 stays like that
    level_t level_      = -1;               //!< the level of the block
    real_t  xyz_[3]     = {0.0, 0.0, 0.0};  //!< the origin of the block
    real_t  hgrid_[3]   = {0.0, 0.0, 0.0};  //!< the grid spacing of the block

    std::map<std::string, mem_ptr> mem_map_;  //<! a map of the pointers to the actual data

    mem_ptr coarse_ptr_;  //!< a memory reserved for coarser version of myself, includes ghost points

    // list of ghosting
    std::list<GhostBlock<GridBlock*>*> local_sibling_;         //<! local neighbors at my resolution
    std::list<GhostBlock<GridBlock*>*> local_parent_;          //!< local neighbors coarser (neighbor to me)
    std::list<GhostBlock<GridBlock*>*> local_parent_reverse_;  //!< local neighbors coarser (me to neighbors)
    std::list<GhostBlock<MPI_Aint>*>   ghost_sibling_;         //<! ghost neighbors at my resolution
    std::list<GhostBlock<MPI_Aint>*>   ghost_parent_;          //!<ghost neighbors coarser (neighbor to me)
    std::list<GhostBlock<MPI_Aint>*>   ghost_children_;        //!<ghost neighbors coarser (neighbor to me)
    std::list<GhostBlock<MPI_Aint>*>   ghost_parent_reverse_;  //!<  ghost neighbors coarser (me to neighbors)
    std::list<PhysBlock*>              phys_;                  //!<  physical boundary condition

    // list of dependency = how to create my information after refinement/coarsening
    sid_t      n_dependency_active_ = 0;
    GridBlock* dependency_[P8EST_CHILDREN];  //!< the pointer to the dependency block

   public:
    GridBlock(const real_t length, const real_t xyz[3], const sid_t level);
    // ~GridBlock();

    /**
     * @name Memory Layout Implementation
     * 
     * @{ */
    inline lid_t gs() const override { return M_GS; }
    inline lid_t stride() const override { return M_STRIDE; }
    inline lid_t start(const int id) const override { return 0; }
    inline lid_t end(const int id) const override { return M_N; }
    /** @} */

    inline sid_t  level() const { return level_; }
    inline real_t xyz(const int id) const { return xyz_[id]; }
    inline real_t hgrid(const int id) const { return hgrid_[id]; }
    const real_t* hgrid() const { return hgrid_; }
    const real_t* xyz() const { return xyz_; }

    /**
     * @name Status level management
     *
     * @{ */
    sid_t status_level() const { return status_lvl_; };
    void  ResetStatus() { status_lvl_ = 0; };
    void  UpdateStatusCriterion(m_ptr<const Wavelet> interp, const real_t rtol, const real_t ctol, m_ptr<const Field> field_citerion, m_ptr<Prof> profiler);
    void  ComputeDetails(m_ptr<const Wavelet> interp, m_ptr<const Field> criterion, m_ptr<const Field> details);
    /** @} */

    /**
     * @name datamap access
     * @{
     */
    data_ptr data(m_ptr<const Field> fid, const lda_t ida = 0) const;
    mem_ptr  pointer(m_ptr<const Field> fid, const lda_t ida = 0) const;
    /** @} */

    /**
     * @name field management
     * 
     * @{
     */
    void AddField(m_ptr<Field> fid);
    void DeleteField(m_ptr<const Field> fid);
    void AddFields(const std::map<std::string, m_ptr<Field> >* fields);
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
    /** @} */

    /**
     * @name coarse pointer management
     * @{
     */
    mem_ptr coarse_ptr() const { return coarse_ptr_; }
    // void    coarse_ptr(mem_ptr ptr) { coarse_ptr_ = ptr; }
    void AllocateCoarsePtr(const size_t memsize);
    /**@} */

    /**
     * @name return the ghost lists
     * @{
     */
    std::list<GhostBlock<GridBlock*>*>* local_sibling() { return &local_sibling_; };
    std::list<GhostBlock<GridBlock*>*>* local_parent() { return &local_parent_; };
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

static inline GridBlock* p4est_GetGridBlock(qdrt_t* quad) {
    return *(reinterpret_cast<GridBlock**>(quad->p.user_data));
}

static inline void p4est_SetGridBlock(qdrt_t* quad, GridBlock* block) {
    *(reinterpret_cast<GridBlock**>(quad->p.user_data)) = block;
}

#endif  // SRC_GRIDBLOCK_HPP_
