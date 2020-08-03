#ifndef SRC_GHOST_HPP_
#define SRC_GHOST_HPP_

#include <list>

#include "doop.hpp"
#include "field.hpp"
#include "forestgrid.hpp"
#include "ghostblock.hpp"
#include "interpolator.hpp"
#include "murphy.hpp"
#include "physblock.hpp"
#include "prof.hpp"

// #define M_NGHOST (M_N * M_N * M_N)

// alias boring names
using GBLocal      = GhostBlock<GridBlock *>;
using GBMirror     = GhostBlock<MPI_Aint>;
using GBPhysic     = PhysBlock;
using ListGBLocal  = std::list<GBLocal *>;
using ListGBMirror = std::list<GBMirror *>;
using listGBPhysic = std::list<GBPhysic *>;

class Ghost;

/**
 * @brief pointer to an member function of the class @ref Ghost
 */
using gop_t = void (Ghost::*)(const qid_t *qid, GridBlock *block,const Field *fid);

/**
 * @brief given an associated grid, performs the ghost update for a given field, in a given component
 * 
 * This class is highly linked to our use of the wavelet, please have a look at @ref doc/ghost.md
 *  
 * @warning the procedure assumes that the grid is correctly balanced!
 * 
 */
class Ghost {
   protected:
    lda_t    ida_           = -1;                  //!< current ghosting dimension
    level_t  min_level_     = -1;                  //!< minimum active level, min_level included
    level_t  max_level_     = P8EST_MAXLEVEL + 1;  //!< maximum active level, max_level included
    iblock_t n_active_quad_ = -1;                  //!< the number of quadrant that needs to have ghost informations

    MPI_Group mirror_origin_group_ = MPI_GROUP_NULL;  //!< group of ranks that will emit/origin a RMA to access my mirrors
    MPI_Group mirror_target_group_ = MPI_GROUP_NULL;  //!< group of ranks that will be targeted by my RMA calls to access mirrors
    iblock_t  n_mirror_to_send_    = 0;               //!< get how many mirrors to send
    MPI_Win   mirrors_window_      = MPI_WIN_NULL;    //!< MPI Window for the RMA communication
    real_t *  mirrors_             = nullptr;         //!< memory space for the mirror blocks, computed using n_mirror_to_send_
    MPI_Win   local2disp_window_   = MPI_WIN_NULL;    //!< MPI Window for the RMA communication (non null only during the initlist function)
    MPI_Aint *local2disp_          = nullptr;         //!< for each quadrant, indicate its corresponding mirror ID (non null only during the initlist function)

    ForestGrid *  grid_;    //!< pointer to the associated @ref ForestGrid, shared, not owned
    Interpolator *interp_;  //!< pointer to the associated @ref Interpolator, shared, not owned

    ListGBLocal ** block_sibling_;         //!<  list of local block on my resolution
    ListGBLocal ** block_parent_;          //!<  list of local block coarser (neighbors to me)
    ListGBLocal ** block_parent_reverse_;  //!<  list of local block coarser (me to neighbors)
    ListGBMirror **ghost_sibling_;         //!<  list of mirror ghosts on my resolution
    ListGBMirror **ghost_children_;        //!<  list of mirror ghosts finer than me
    ListGBMirror **ghost_parent_;          //!<  list of mirror ghost coarser (neighbors to me)
    ListGBMirror **ghost_parent_reverse_;  //!<  list of mirror ghost coarser (neighbors to me)
    listGBPhysic **phys_;                  //!<  physical boundary condition

   public:
    Ghost(ForestGrid *grid, Interpolator *interp);
    Ghost(ForestGrid *grid, const level_t min_level, const level_t max_level, Interpolator *interp);
    ~Ghost();

    /**
     * @name RMA-based high-level ghosting - post and wait
     * @{
     */
    void PullGhost_Post(Field *field, const lda_t ida);
    void PullGhost_Wait(Field *field, const lda_t ida);
    /** @}*/

    /**
     * @name RMA-based low-level ghosting - get and put the values
     * @{
     */
    void InitList4Block(const qid_t *qid,       GridBlock *block);
    void PushToWindow4Block(const qid_t *qid,   GridBlock *block, const Field *fid);
    void GetGhost4Block_Post(const qid_t *qid,  GridBlock *block, const Field *fid);
    void GetGhost4Block_Wait(const qid_t *qid,  GridBlock *block, const Field *fid);
    void PutGhost4Block_Post(const qid_t *qid,  GridBlock *block, const Field *fid);
    void PutGhost4Block_Wait(const qid_t *qid,  GridBlock *block, const Field *fid);
    void PullFromWindow4Block(const qid_t *qid, GridBlock *block, const Field *fid);
    /** @}*/

   protected:

    /**
     * @name Init and free the comm/list
     * @{
     */
    void InitComm_();
    void FreeComm_();
    void InitList_();
    void FreeList_();
    /** @}*/

    /**
     * @name Different operations needed on the block itself, given a ghost list
     * @{
     */
    inline void Compute4Block_Myself2Coarse_(const qid_t *qid, GridBlock *cur_block, const Field *fid, mem_ptr ptr_trg);
    inline void Compute4Block_Copy2Myself_(const ListGBLocal *ghost_list, const Field *fid, data_ptr data_trg);
    inline void Compute4Block_Copy2Coarse_(const ListGBLocal *ghost_list, const Field *fid, mem_ptr ptr_trg);
    inline void Compute4Block_GetRma2Myself_(const ListGBMirror *ghost_list, const Field *fid, data_ptr data_trg);
    inline void Compute4Block_GetRma2Coarse_(const ListGBMirror *ghost_list, const Field *fid, mem_ptr ptr_trg);
    inline void Compute4Block_Refine_(const ListGBLocal *ghost_list, const mem_ptr ptr_src, data_ptr data_trg);
    inline void Compute4Block_Refine_(const ListGBMirror *ghost_list, const mem_ptr ptr_src, data_ptr data_trg);
    inline void Compute4Block_Coarsen2Coarse_(data_ptr data_src, data_ptr ptr_trg);
    inline void Compute4Block_Copy2Parent_(const ListGBLocal *ghost_list, const mem_ptr ptr_src, const Field *fid);
    inline void Compute4Block_PutRma2Parent_(const ListGBMirror *ghost_list, const mem_ptr ptr_src);
    inline void Compute4Block_Phys2Myself_(const qid_t *qid, GridBlock *cur_block, const Field *fid);
    /** @}*/

    void LoopOnMirrorBlock_(const gop_t op, const Field *field);
};

#endif  // SRC_GHOST_HPP_
