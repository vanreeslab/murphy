#ifndef __GHOST_HPP
#define __GHOST_HPP

#include <p8est.h>
#include <p8est_mesh.h>

#include "ghostblock.hpp"
#include "grid.hpp"
#include "murphy.hpp"
#include "operator.hpp"
#include "physblock.hpp"

using std::list;

#define M_NGHOST (M_N * M_N * M_N)
#define M_NNEIGHBOR 26

/**
 * @brief the Ghost class is both a constant operator (see @ref ConstOperatorF) when performing the ghost exchange
 * and a simple Operator (see @ref OperatorS) when performing the initialization.
 * 
 * 
 */
class Ghost : public ConstOperatorF, public OperatorS {
   protected:
    list<GhostBlock*>** block_sibling_; //!<  list of blocks that are finer or same resolution 
    list<GhostBlock*>** block_parent_;  //!<  list of blocks that are coarser
    list<GhostBlock*>** ghost_sibling_; //!<  list of ghosts that are finer or same resolution 
    list<GhostBlock*>** ghost_parent_;  //!<  list of ghosts that are coarser 
    list<PhysBlock*>**  phys_;          //!<  physical blocks 

    real_t* mirrors_ = nullptr;
    real_t* ghosts_  = nullptr;

    Grid* grid_ = nullptr;

   public:
    Ghost(Grid* grid);
    ~Ghost();

    void reset();

    /**
     *  @defgroup Ghosts initialization
     * 
     * Implement the function of the @ref OperatorS class 
     *  @{
     */
    void apply(const qid_t* qid, GridBlock* block) override;
    // void operator()(Grid* grid) override;
    /** @} */

    /**
     *  @defgroup Ghosts application 
     * 
     * Implement the function of the @ref ConstOperatorF class 
     * 
     *  @{
     */
    void apply(const qid_t* qid, GridBlock* block, const Field* fid) override;
    // void operator()(Grid* grid, Field* field) override;
    /** @} */

   protected:
    //
};

#endif