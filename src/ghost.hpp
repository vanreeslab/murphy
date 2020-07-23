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
    iblock_t n_active_quad_ = -1;                   //!< the number of quadrant that needs to have ghost informations
    sid_t    nghost_[2]     = {0, 0};              //!< the number of ghost (front,back) that are actually needed

    // information that tracks which block is involved
    iblock_t  n_mirror_to_send_ = 0;        //!< get how many mirrors to send
    iblock_t  n_ghost_to_recv_  = 0;        //!< get how many ghosts to recv
    iblock_t *local_to_mirrors  = nullptr;  //!< for each registered mirror to send, stores the p4est mirror id
    iblock_t *ghost_to_local_   = nullptr;  //!< for each p4est ghost block stores the local id of the ghost to recv

    lid_t        n_send_request_ = 0;        //!< the number of send requests by level
    lid_t        n_recv_request_ = 0;        //!< the number of receive requests by level
    MPI_Request *mirror_send_    = nullptr;  //!< the send requests for the mirrors
    MPI_Request *ghost_recv_     = nullptr;  //!< the receive requests for the ghosts

    //---------------
    // RMA
    MPI_Group mirror_origin_group_ = MPI_GROUP_NULL;  //!< group of ranks that will emit a RMA to access my mirrors
    MPI_Group mirror_target_group_ = MPI_GROUP_NULL;  //!< group of ranks that will be target by my RMA calls to access mirrors

    // acess to the mirror data
    MPI_Win mirrors_window_;     //!< MPI Window for the RMA communication
    real_t *mirrors_ = nullptr;  //!< memory space for the mirror blocks, computed using n_mirror_to_send_

    // access to the mirror displacement -> non null only during the initlist function
    MPI_Win   local2disp_window_;     //!< MPI Window for the RMA communication
    MPI_Aint *local2disp_ = nullptr;  //!< for each quadrant, indicate its corresponding mirror ID

    //---------------

    real_t *ghosts_         = nullptr;       //!< memory space for the ghost blocks

    ForestGrid *  grid_;        //!< pointer to the associated @ref ForestGrid, shared, not owned
    real_p *      coarse_tmp_;  //!< working memory that contains a coarse version of the current block, one per thread
    Interpolator *interp_;      //!< pointer to the associated @ref Interpolator, shared, not owned

    ListGBLocal ** block_sibling_;   //!<  list of blocks that are on the  same resolution
    ListGBLocal ** block_children_;  //!<  list of blocks that are finer
    ListGBLocal ** block_parent_;    //!<  list of blocks that are coarser
    ListGBMirror **ghost_sibling_;   //!<  list of ghosts that are on the same resolution
    ListGBMirror **ghost_children_;  //!<  list of ghosts that are on the same resolution
    ListGBMirror **ghost_parent_;    //!<  list of ghosts that are coarser
    listGBPhysic **phys_;            //!<  physical blocks

   public:
    Ghost(ForestGrid *grid, Interpolator *interp);
    Ghost(ForestGrid *grid, const level_t min_level, const level_t max_level, Interpolator *interp);
    ~Ghost();

    /**
     * @name RMA-based ghosting -- divided in 4 steps (order matters!)
     * @{
     */
    void GetGhost_Post(Field *field, const sid_t ida);  // step 1
    void GetGhost_Wait(Field *field, const sid_t ida);  // step 2
    void PutGhost_Post(Field *field, const sid_t ida);  // step 3
    void PutGhost_Wait(Field *field, const sid_t ida);  // step 4
    // block functions
    void GetGhost4Block_Post(const qid_t *qid, GridBlock *block, Field *fid);  // step 1
    void GetGhost4Block_Wait(const qid_t *qid, GridBlock *block, Field *fid);  // step 2
    void PutGhost4Block_Post(const qid_t *qid, GridBlock *block, Field *fid);  // step 3
    void PutGhost4Block_Wait(const qid_t *qid, GridBlock *block, Field *fid);  // step 4
    /** @}*/

    void PushToMirror(Field *field, const sid_t ida);
    void MirrorToGhostSend(Prof *prof);
    void MirrorToGhostRecv(Prof *prof);
    void PullFromGhost(Field *field, const sid_t ida);

    /**
     *  @name Execute on each block
     * 
     *  @{
     */
    void InitList4Block(const qid_t *qid, GridBlock *block);
    void PushToWindow4Block(const qid_t *qid, GridBlock *block, Field *fid);

    void PushToMirror4Block(const qid_t *qid, GridBlock *block, Field *fid);
    void PullFromGhost4Block(const qid_t *qid, GridBlock *cur_block, Field *fid);
    /** @} */

   protected:
    // void InitComm_();
    void InitComm_();
    void FreeComm_();
    void InitList_();
    void FreeList_();

    void Compute4Block_Copy2Myself_(const ListGBLocal *ghost_list, Field *fid, GridBlock *block_trg, real_t *data_trg);
    void Compute4Block_Copy2Coarse_(const ListGBLocal *ghost_list, Field *fid, GridBlock *block_trg, real_t *ptr_trg);
    void Compute4Block_GetRma2Myself_(const ListGBMirror *ghost_list, Field *fid, GridBlock *block_trg, real_t *data_trg);
    void Compute4Block_GetRma2Coarse_(const ListGBMirror *ghost_list, Field *fid, GridBlock *block_trg, real_t *ptr_trg);
    void Compute4Block_Myself2Coarse_(const qid_t *qid, GridBlock *cur_block, Field *fid, real_t *ptr_trg);
    void Compute4Block_Refine_(const ListGBLocal *ghost_list, real_t *ptr_src, real_t *data_trg);
    void Compute4Block_Refine_(const ListGBMirror *ghost_list, real_t *ptr_src, real_t *data_trg);
    void Compute4Block_Phys2Myself_(const qid_t *qid, GridBlock *cur_block, Field *fid);

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
