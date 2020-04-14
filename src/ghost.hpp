#ifndef __GHOST_HPP
#define __GHOST_HPP

#include <p8est.h>
#include <p8est_mesh.h>

#include "ghostblock.hpp"
#include "grid.hpp"
#include "murphy.hpp"
#include "operator.hpp"
#include "physblock.hpp"
#include "interpolator.hpp"

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

    Grid* grid_ = nullptr;  //!< the associated grid

    sid_t ida_ = -1;  //!< current ghosting dimension

    real_p* coarse_tmp_; //!< working memory that contains a coarse version of the current block, one per thread

    Interpolator* interpolator_;//!< current interpolator used

   public:
    Ghost(Grid* grid,Interpolator* interpolator);
    ~Ghost();

    /**
     * @brief apply the ghost computation for a field in a given direction
     * 
     * @param field the current field
     * @param current_ida the working dimension
     */
    void pull(Field* field);

    /**
     *  @name OperatorS implementation
     *  @{
     */
    void apply(const qid_t* qid, GridBlock* block) override;
    /** @} */

    /**
     *  @name ConstOperatorF implementation
     *  @{
     */
    void apply(const qid_t* qid, GridBlock* block, Field* fid) override;
    /** @} */

   protected:
    //
};

#endif