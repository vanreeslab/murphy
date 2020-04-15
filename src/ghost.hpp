#ifndef SRC_GHOST_HPP_
#define SRC_GHOST_HPP_

#include <list>

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

    ForestGrid* grid_;
    Interpolator* interp_;

    real_p* coarse_tmp_;  //!< working memory that contains a coarse version of the current block, one per thread

   public:
    Ghost(ForestGrid* grid);
    ~Ghost();

    
    void Pull(Field* field,Interpolator* interp);

    /**
     *  @name OperatorS implementation
     *  @{
     */
    void ApplyOperatorS(const qid_t* qid, GridBlock* block) override;
    /** @} */

    /**
     *  @name ConstOperatorF implementation
     *  @{
     */
    void ApplyOperatorF(const qid_t* qid, GridBlock* block, Field* fid) override;
    /** @} */

   protected:
    //
};

#endif  // SRC_GHOST_HPP_