#ifndef SRC_GHOST_HPP_
#define SRC_GHOST_HPP_

#include <list>

#include "field.hpp"
#include "forestgrid.hpp"
#include "ghostblock.hpp"
#include "interpolator.hpp"
#include "murphy.hpp"
#include "operator.hpp"
#include "physblock.hpp"
#include "prof.hpp"

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
 * It is both a constant operator (see @ref OperatorF) when performing the ghost exchange
 * and a simple Operator (see @ref OperatorS) when performing the initialization.
 * 
 * Ghosting relies on two structures, the blocks (local leafs) and the ghosts (other rank's leafs).
 * Only the non-local ghosts will have to be received. In the origin tree, they are called mirrors.
 * It means that a mirror is a ghost for one or multiple rank and a ghost is the block once received.
 * Then we interpolate (or copy) the information to the correct spot.
 * 
 * For the moment, all the block is being send/received. One optimisation is to only send a layer of 2x4 points in each direction.
 * Indeed, because we might have a coarser neighbor, we need to coarsen twice a block (worst case) and then need 4 points in every direction
 * 
 */
class Ghost : public OperatorF, public OperatorS {
   protected:
    list<GhostBlock*>** block_sibling_;  //!<  list of blocks that are finer or same resolution
    list<GhostBlock*>** block_parent_;   //!<  list of blocks that are coarser
    list<GhostBlock*>** ghost_sibling_;  //!<  list of ghosts that are finer or same resolution
    list<GhostBlock*>** ghost_parent_;   //!<  list of ghosts that are coarser
    list<PhysBlock*>**  phys_;           //!<  physical blocks

    lid_t   n_mirror_to_send_ = 0;        //!< get how many mirrors we have to be send, as a sum of the mirrors send for each rank (one mirror can be send to multiple ranks!!)
    lid_t*  mirrors_to_local_ = nullptr;  //!< defines the relation between a send-mirror and it's local ID used to access p4est mirrors array (one p4est mirrors quad can be send multiple times to different ranks)
    real_t* mirrors_          = nullptr;  //!< memory space for the mirror blocks, computed using n_mirror_to_send_
    real_t* ghosts_           = nullptr;  //!< memory space for the ghost blocks

    sid_t ida_ = -1;  //!< current ghosting dimension

    lid_t        n_send_request_ = 0;        //!< the number of send requests
    lid_t        n_recv_request_ = 0;        //!< the number of receive requests
    MPI_Request* mirror_send_    = nullptr;  //!< the send requests for the mirrors
    MPI_Request* ghost_recv_     = nullptr;  //!< the receive requests for the ghosts

    ForestGrid*   grid_;    //!< pointer to the associated @ref ForestGrid, shared, not owned
    Interpolator* interp_;  //!< pointer to the associated @ref Interpolator, shared, not owned

    real_p* coarse_tmp_;  //!< working memory that contains a coarse version of the current block, one per thread

   public:
    Ghost(ForestGrid* grid, Interpolator* interp);
    ~Ghost();

    void PushToMirror(Field* field, sid_t ida);
    void MirrorToGhostSend(Prof* prof);
    void MirrorToGhostRecv(Prof* prof);
    void PullFromGhost(Field* field, sid_t ida);

    /**
     *  @name OperatorS implementation
     *  @{
     */
    void ApplyOpS(const qid_t* qid, GridBlock* block) override;
    /** @} */

    /**
     *  @name ConstOperatorF implementation
     * 
     *  @{
     */
    void ApplyOpF(const qid_t* qid, GridBlock* block, Field* fid) override;
    /** @} */

   protected:
    void InitComm_();
    void InitList_(const qid_t* qid, GridBlock* block);

    void PushToMirror_(const qid_t* qid, GridBlock* block, Field* fid);
    void PullFromGhost_(const qid_t* qid, GridBlock* block, Field* fid);

    void LoopOnMirrorBlock_(const gop_t op, Field* field);
};

/**
 * @brief Given a face, edge or a corner returns the outgoing normal (sign)
 * 
 * @param ibidule for a face (`ibidule<6`), an edge (`6<= ibidule < 18`) or a corner (`18<= ibidule < 26`)
 * @param sign the sign of the outgoing normal
 */
static void GhostGetSign(sid_t ibidule, real_t sign[3]) {
    // we need to find the sign = the direction of the normal:
    sign[0] = 0;
    sign[1] = 0;
    sign[2] = 0;

    // check depending on the plane, the edge of the corner
    if (ibidule < 6) {
        sid_t dir = ibidule / 2;
        sign[dir] = ((ibidule % 2) == 1) ? 1 : -1;
    } else if (ibidule < 18) {
        sid_t iedge = ibidule - 6;
        /*
        the plane convention for the sign variable convention for the sign
        2 +--------------+ 3
          |              |
          |              |

          |dir2          |
          |              |
        0 +--------------+ 1
            dir1
        */
        sid_t dir  = iedge / 4;           // this is the direction of the edge
        sid_t dir1 = (dir == 0) ? 1 : 0;  // dir1 in the plane: dir1 = x if dir = y or z or y if dir = x
        sid_t dir2 = (dir == 2) ? 1 : 2;  // dir2 in the plane: dir2 = y if dir=z, = z if dir=x or dir = y
        sign[dir1] = ((iedge % 4) % 2) == 1 ? +1 : -1;
        sign[dir2] = ((iedge % 4) / 2) == 1 ? +1 : -1;
    } else {
        sid_t icorner = ibidule - 18;
        sign[0]       = (icorner % 2) == 1 ? +1 : -1;
        sign[1]       = ((icorner % 4) / 2) == 1 ? +1 : -1;
        sign[2]       = (icorner / 4) == 1 ? +1 : -1;
    }
}

/**
 * @brief given a starting id in a block and a number of coarse ghost points, returns the corresponding starting id for the coarse representation of it
 * 
 * @note this is magic...
 * 
 * we compute:
 * b = a + M_N, such that
 *      if a is in the ghost points, b < M_N
 *      if a is in the center points, M_N <= b < 2* M_N
 *      if a is in the ghost points, 2* M_N <= b
 * 
 * c is binary: 0 if a is in the ghost points and 1 if a is in the center points
 * 
 * d is the number of points taken into the ghost points or the number of points interior to the block
 * 
 * e is the scaled d to the coarse block:
 *  if d is interior: c=1, and we divide by 2 the id
 *  if d is the number in the ghost points:  c=0 and
 *      if d contains all the ghost points, we scale it to contain all the coarse gp
 *      if d contains a part of the ghost points, we abort
 * 
 * we return a sum of:
 * (a / M_N) * M_HN = the rescaled length of a to the new length (M_HN):
 *      if a < M_N, it will be 0
 *      if a >= M_N, it will be 1
 * (a % M_N) / (c + 1) = the rest of the lentgh, divided by c+1=2 if a is in the center or =1 if a is on the side
 * 
 * As a summary:
 *      if a is in the ghost points, a > 0, c = 0 =>  returns 0 + a/1
 *      if a is in the center points, 0 <= a < M_N => returns 0 + a/2
 *      if a is in the ghost points, a >= M_N => returns M_HN + (a%M_N)/1 
 * 
 * @param a the id in the block
 * @param cgs the number of ghost points in the coarse level
 * @return lid_t 
 */
static inline lid_t CoarseFromBlock(const lid_t a, const lid_t cgs) {
    const lid_t b = (a + M_N);
    const lid_t c = (b / M_N) % 2;
    const lid_t d = (a % M_N);
    const lid_t e = c * (d / 2) + (1 - c) * ((d / M_GS) * cgs);
    m_assert(((1 - c) * (d % M_GS)) == 0, "this should NOT happen");
    return (a / M_N) * M_HN + e;
}

#endif  // SRC_GHOST_HPP_