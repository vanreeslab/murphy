#ifndef SRC_GRIDBLOCK_HPP_
#define SRC_GRIDBLOCK_HPP_

#include <mpi.h>

#include <limits>
#include <list>
#include <string>
#include <unordered_map>

#include "defs.hpp"
#include "field.hpp"
#include "ghostblock.hpp"
#include "interpolator.hpp"
#include "memlayout.hpp"
#include "murphy.hpp"
#include "p8est.h"
#include "physblock.hpp"
#include "wavelet.hpp"
#include "forestgrid.hpp"


/**
 * @brief implements a @ref Block that is used as a leaf for the tree
 * 
 */
class GridBlock : public MemLayout {
   protected:
    bool    lock_     = false;            //!< lock the block, indicating that no refinement/coarsening can happen on the block
    level_t level_    = -1;               //!< the level of the block
    real_t  xyz_[3]   = {0.0, 0.0, 0.0};  //!< the origin of the block
    real_t  hgrid_[3] = {0.0, 0.0, 0.0};  //!< the grid spacing of the block

    std::unordered_map<std::string, mem_ptr> data_map_;  //<! a map of the pointers to the actual data

    mem_ptr coarse_ptr_ = nullptr;  //!< a memory reserved for coarser version of myself, includes ghost points

    // list of ghosting
    std::list<GhostBlock<GridBlock*>*> local_sibling_;         //<! local neighbors at my resolution
    std::list<GhostBlock<GridBlock*>*> local_parent_;          //!< local neighbors coarser (neighbor to me)
    std::list<GhostBlock<GridBlock*>*> local_parent_reverse_;  //!< local neighbors coarser (me to neighbors)
    std::list<GhostBlock<MPI_Aint>*>   ghost_sibling_;         //<! ghost neighbors at my resolution
    std::list<GhostBlock<MPI_Aint>*>   ghost_parent_;          //!<ghost neighbors coarser (neighbor to me)
    std::list<GhostBlock<MPI_Aint>*>   ghost_children_;        //!<ghost neighbors coarser (neighbor to me)
    std::list<GhostBlock<MPI_Aint>*>   ghost_parent_reverse_;  //!<  ghost neighbors coarser (me to neighbors)
    std::list<PhysBlock*>              phys_;                  //!<  physical boundary condition

   public:
    GridBlock(const real_t length, const real_t xyz[3], const sid_t level);
    ~GridBlock();

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
     * @name Lock management
     * 
     * @{ */
    void lock() { lock_ = true; }
    void unlock() { lock_ = false; }
    bool locked() const { return lock_; }
    /**@} */

    /**
     * @name datamap access
     * 
     * @{
     */
    // data = memory address of (0,0,0)
    data_ptr data(const Field* fid);
    data_ptr data(const Field* fid, const sid_t ida);
    // pointer = raw data pointe
    mem_ptr pointer(const Field* fid);
    mem_ptr pointer(const Field* fid, const sid_t ida);
    /** @} */

    /**
     * @name field management
     * 
     * @{
     */
    void AddField(Field* fid);
    void DeleteField(Field* fid);
    void AddFields(const std::unordered_map<std::string, Field*>* fields);
    /** @} */

    /**
     * @name handle the ghost data pointer
     * @{
     */
    mem_ptr coarse_ptr() const { return coarse_ptr_; }
    void    coarse_ptr(mem_ptr ptr) { coarse_ptr_ = ptr; }
    void    AllocateCoarsePtr(const size_t memsize);
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

    void GhostInitLists(const qid_t* qid, const ForestGrid* grid, const InterpolatingWavelet* interp, MPI_Win local2disp_window);
    void GhostFreeLists();
    void GhostGet_Post(const Field* field, const lda_t ida, const InterpolatingWavelet* interp, MPI_Win mirrors_window);
    void GhostGet_Wait(const Field* field, const lda_t ida, const InterpolatingWavelet* interp);
    void GhostPut_Post(const Field* field, const lda_t ida, const InterpolatingWavelet* interp, MPI_Win mirrors_window);
    void GhostPut_Wait(const Field* field, const lda_t ida, const InterpolatingWavelet* interp);
};

#endif  // SRC_GRIDBLOCK_HPP_
