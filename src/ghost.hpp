#ifndef SRC_GHOST_HPP_
#define SRC_GHOST_HPP_

#include <list>

#include "field.hpp"
#include "forestgrid.hpp"
#include "ghostblock.hpp"
#include "interpolator.hpp"
#include "murphy.hpp"
#include "doop.hpp"
#include "physblock.hpp"
#include "prof.hpp"

using std::list;

// #define M_NGHOST (M_N * M_N * M_N)

// alias boring names
using GBLocal      = GhostBlock<GridBlock *>;
using GBMirror     = GhostBlock<MPI_Aint>;
using GBPhysic     = PhysBlock;
using ListGBLocal  = list<GBLocal *>;
using ListGBMirror = list<GBMirror *>;
using listGBPhysic = list<GBPhysic *>;

class Ghost;

/**
 * @brief pointer to an member function of the class @ref Ghost
 */
using gop_t = void (Ghost::*)(const qid_t *qid, GridBlock *block, Field *fid);

/**
 * @brief performs the ghost update of a given grid field, in one dimension
 * 
 * For some description of the p4est library related to the ghost management, see @ref doc/p4est.md
 * 
 * Ghosting relies on two structures, the blocks (local leafs) and the ghosts (other rank's leafs).
 * Only the non-local ghosts will have to be received. In the origin tree, they are called mirrors.
 * It means that a mirror is a ghost for one or multiple rank and a ghost is the block once received.
 * Then we interpolate (or copy) the information to the correct spot.
 * 
 * For the moment, all the block is being send/received. One optimisation is to only send an inner layer of points in each direction.
 * 
 * @warning the procedure assumes that the grid is correctly balanced!
 * 
 */

class Ghost {
   protected:
    sid_t    ida_           = -1;                  //!< current ghosting dimension
    level_t  min_level_     = -1;                  //!< minimum active level, min_level included
    level_t  max_level_     = P8EST_MAXLEVEL + 1;  //!< maximum active level, max_level included
    iblock_t n_active_quad_ = -1;                  //!< the number of quadrant that needs to have ghost informations
    sid_t    nghost_[2]     = {0, 0};              //!< the number of ghost (front,back) that are actually needed

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
    void PullGhost_Post(Field *field, const sid_t ida);
    void PullGhost_Wait(Field *field, const sid_t ida);
    /** @}*/

    /**
     * @name RMA-based low-level ghosting - get and put the values
     * @{
     */
    void InitList4Block(const qid_t *qid, GridBlock *block);
    void PushToWindow4Block(const qid_t *qid, GridBlock *block, Field *fid);
    void GetGhost4Block_Post(const qid_t *qid, GridBlock *block, Field *fid);
    void GetGhost4Block_Wait(const qid_t *qid, GridBlock *block, Field *fid);
    void PutGhost4Block_Post(const qid_t *qid, GridBlock *block, Field *fid);
    void PutGhost4Block_Wait(const qid_t *qid, GridBlock *block, Field *fid);
    void PullFromWindow4Block(const qid_t *qid, GridBlock *block, Field *fid);
    /** @}*/

   protected:
    // void InitComm_();
    void InitComm_();
    void FreeComm_();
    void InitList_();
    void FreeList_();

    inline void Compute4Block_Myself2Coarse_(const qid_t *qid, GridBlock *cur_block, Field *fid, real_t *ptr_trg);
    inline void Compute4Block_Copy2Myself_(const ListGBLocal *ghost_list, Field *fid, GridBlock *block_trg, real_t *data_trg);
    inline void Compute4Block_Copy2Coarse_(const ListGBLocal *ghost_list, Field *fid, GridBlock *block_trg, real_t *ptr_trg);
    inline void Compute4Block_GetRma2Myself_(const ListGBMirror *ghost_list, Field *fid, GridBlock *block_trg, real_t *data_trg);
    inline void Compute4Block_GetRma2Coarse_(const ListGBMirror *ghost_list, Field *fid, GridBlock *block_trg, real_t *ptr_trg);
    inline void Compute4Block_Refine_(const ListGBLocal *ghost_list, real_t *ptr_src, real_t *data_trg);
    inline void Compute4Block_Refine_(const ListGBMirror *ghost_list, real_t *ptr_src, real_t *data_trg);
    inline void Compute4Block_Coarsen2Coarse_(real_t *data_src, real_t *ptr_trg);
    inline void Compute4Block_Copy2Parent_(const ListGBLocal *ghost_list, real_t *ptr_src, Field *fid);
    inline void Compute4Block_PutRma2Parent_(const ListGBMirror *ghost_list, real_t *ptr_src, Field *fid);
    inline void Compute4Block_Phys2Myself_(const qid_t *qid, GridBlock *cur_block, Field *fid);

    // void Compute4Block_Sibling_(const list<GhostBlock *> *ghost_list, const bool do_coarse, InterpFunction *copy, GridBlock *cur_block, Field *fid, real_t *coarse_mem);
    // void Compute4Block_FromParent_(const list<GhostBlock *> *ghost_list, InterpFunction *copy, GridBlock *cur_block, Field *fid, real_t *coarse_mem);

    void LoopOnMirrorBlock_(const gop_t op, Field *field);
    void LoopOnGhostBlock_(const gop_t op, Field *field);

    // void PullFromGhost4Block_Children(const qid_t *qid, GridBlock *cur_block, Field *fid,
    //                                  const bool do_coarse, SubBlock *ghost_block, SubBlock *coarse_block, real_t *coarse_mem);
    // void PullFromGhost4Block_Sibling_(const qid_t *qid, GridBlock *cur_block, Field *fid,
    //                                   const bool do_coarse, SubBlock *ghost_block, SubBlock *coarse_block, real_t *coarse_mem);
    void PullFromGhost4Block_FromParent_(const qid_t *qid, GridBlock *cur_block, Field *fid,
                                        const bool do_coarse, SubBlock *ghost_block, SubBlock *coarse_block, real_t *coarse_mem);
    void PullFromGhost4Block_Myself_(const qid_t *qid, GridBlock *cur_block, Field *fid, SubBlock *coarse_block, real_t *coarse_mem);
    void PullFromGhost4Block_ToParent_(const qid_t *qid, GridBlock *cur_block, Field *fid,
                                      SubBlock *ghost_block, SubBlock *coarse_block, real_t *coarse_mem);
    // void PullFromGhost4Block_Physics(const qid_t *qid, GridBlock *cur_block, PhysBlock *gblock, Field *fid, real_t hgrid[3], real_t *data);

    // void CreateOnLevels_(ForestGrid *grid, Interpolator *interp, const level_t min_level, const level_t max_level);
};

#endif  // SRC_GHOST_HPP_
