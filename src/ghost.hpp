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

using std::list;

#define M_NGHOST (M_N * M_N * M_N)
#define M_NNEIGHBOR 26
#define M_CLEN (2 * M_GS + M_HN)

class Ghost;
/**
 * @brief pointer to an member function of the class @ref Ghost
 */
using gop_t = void (Ghost::*)(const qid_t* qid, GridBlock* block, Field* fid);

/**
 * @brief performs the ghost update of a given grid field
 * 
 * It is both a constant operator (see @ref OperatorF) when performing the ghost exchange
 * and a simple Operator (see @ref OperatorS) when performing the initialization.
 * 
 * 
 */
class Ghost : public OperatorF, public OperatorS {
   protected:
    list<GhostBlock*>** block_sibling_;  //!<  list of blocks that are finer or same resolution
    list<GhostBlock*>** block_parent_;   //!<  list of blocks that are coarser
    list<GhostBlock*>** ghost_sibling_;  //!<  list of ghosts that are finer or same resolution
    list<GhostBlock*>** ghost_parent_;   //!<  list of ghosts that are coarser
    list<PhysBlock*>**  phys_;           //!<  physical blocks

    real_t* mirrors_ = nullptr;  //!< memory space for the mirror blocks
    real_t* ghosts_  = nullptr;  //!< memory space for the ghost blocks

    sid_t ida_ = -1;  //!< current ghosting dimension

    lid_t        n_send_request_ = 0;
    lid_t        n_recv_request_ = 0;
    MPI_Request* mirror_send_    = nullptr;
    MPI_Request* ghost_recv_     = nullptr;

    ForestGrid*   grid_;
    Interpolator* interp_;

    real_p* coarse_tmp_;  //!< working memory that contains a coarse version of the current block, one per thread

   public:
    Ghost(ForestGrid* grid);
    ~Ghost();

    void PushToMirror(Field* field, sid_t ida);
    void MirrorToGhostSend();
    void MirrorToGhostRecv();
    void PullFromGhost(Field* field, sid_t ida, Interpolator* interp);

    /**
     *  @name OperatorS implementation
     *  @{
     */
    void ApplyOperatorS(const qid_t* qid, GridBlock* block) override;
    /** @} */

    /**
     *  @name ConstOperatorF implementation
     * 
     *  @{
     */
    void ApplyOperatorF(const qid_t* qid, GridBlock* block, Field* fid) override;
    /** @} */

   protected:
    void InitComm_();
    void InitList_(const qid_t* qid, GridBlock* block);

    void PushToMirror_(const qid_t* qid, GridBlock* block, Field* fid);
    void PullFromGhost_(const qid_t* qid, GridBlock* block, Field* fid);

    void LoopOnMirrorBlock_(const gop_t op, Field* field);
};



/**
 * @brief given a starting id in a block, returns the corresponding starting id for the coarse representation of it
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
 * we return a sum of:
 * (a / M_N) * M_HN = the rescaled length of a to the new length (M_HN):
 *      if a < M_N, it will be 0
 *      if a >= M_N, it will be 1
 * (a % M_N) / (c + 1) = the rest of the lentgh, divided by 2 if a is in the center.
 * 
 * As a summary:
 *      if a is in the ghost points, a > 0, c = 0 =>  returns 0 + a/1
 *      if a is in the center points, 0 <= a < M_N => returns 0 + a/2
 *      if a is in the ghost points, a >= M_N => returns M_HN + (a%M_N)/1 
 * 
 * @param a 
 * @return lid_t 
 */
static inline lid_t CoarseFromBlock(const lid_t a) {
    const lid_t b = (a + M_N);
    const lid_t c = ( b/ M_N) % 2;
    return (a / M_N) * M_HN + (a % M_N) / (c + 1);
};

#endif  // SRC_GHOST_HPP_