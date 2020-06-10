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
#include "mempool.hpp"

using std::list;

#define M_NGHOST (M_N * M_N * M_N)
#define M_NNEIGHBOR 26

class Ghost;

/**
 * @brief pointer to an member function of the class @ref Ghost
 */
using gop_t = void (Ghost::*)(const qid_t* qid, GridBlock* block, Field* fid);

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
    sid_t    ida_          = -1;                  //!< current ghosting dimension
    level_t  min_level_    = -1;                  //!< minimum active level, min_level included
    level_t  max_level_    = P8EST_MAXLEVEL + 1;  //!< maximum active level, max_level included
    iblock_t n_active_quad_ = 0;                   //!< the number of quadrant that needs to have ghost informations

    // information that tracks which block is involved
    iblock_t  n_mirror_to_send_ = 0;        //!< get how many mirrors to send
    iblock_t  n_ghost_to_recv_  = 0;        //!< get how many ghosts to recv
    iblock_t* local_to_mirrors = nullptr;  //!< for each registered mirror to send, stores the p4est mirror id
    iblock_t* ghost_to_local_   = nullptr;  //!< for each p4est ghost block stores the local id of the ghost to recv

    lid_t        n_send_request_   = 0;  //!< the number of send requests by level
    lid_t        n_recv_request_   = 0;  //!< the number of receive requests by level
    MPI_Request* mirror_send_      = nullptr;  //!< the send requests for the mirrors
    MPI_Request *ghost_recv_       = nullptr;  //!< the receive requests for the ghosts

    real_t *mirrors_ = nullptr;  //!< memory space for the mirror blocks, computed using n_mirror_to_send_
    real_t *ghosts_  = nullptr;  //!< memory space for the ghost blocks

    ForestGrid *  grid_;    //!< pointer to the associated @ref ForestGrid, shared, not owned
    
    // real_p *coarse_tmp_;  //!< working memory that contains a coarse version of the current block, one per thread
    MemPool* mem_pool_; //!< the current memory pool to use for the ghost, not owned
    Interpolator *interp_;  //!< pointer to the associated @ref Interpolator, shared, not owned


    list<GhostBlock *> **block_sibling_;  //!<  list of blocks that are finer or same resolution
    list<GhostBlock *> **block_parent_;   //!<  list of blocks that are coarser
    list<GhostBlock *> **ghost_sibling_;  //!<  list of ghosts that are finer or same resolution
    list<GhostBlock *> **ghost_parent_;   //!<  list of ghosts that are coarser
    list<PhysBlock *> ** phys_;           //!<  physical blocks

   public:
    Ghost(ForestGrid *grid, const level_t min_level, const level_t max_level);
    Ghost(ForestGrid *grid) : Ghost(grid, -1, P8EST_MAXLEVEL + 1){};
    ~Ghost();

    void PushToMirror(Field *field, const sid_t ida);
    void MirrorToGhostSend(Prof *prof);
    void MirrorToGhostRecv(Prof *prof);
    void PullFromGhost(Field* field, const sid_t ida, Interpolator* interp, MemPool* mem_pool);

    /**
     *  @name Execute on each block
     * 
     *  @{
     */
    void InitList4Block(const qid_t *qid, GridBlock *block);
    void PushToMirror4Block(const qid_t *qid, GridBlock *block, Field *fid);
    void PullFromGhost4Block(const qid_t *qid, GridBlock *block, Field *fid);
    /** @} */

   protected:
    void InitComm_();
    void LoopOnMirrorBlock_(const gop_t op, Field *field);
    void LoopOnGhostBlock_(const gop_t op, Field *field);
    // void CreateOnLevels_(ForestGrid *grid, Interpolator *interp, const level_t min_level, const level_t max_level);
};

#endif  // SRC_GHOST_HPP_
