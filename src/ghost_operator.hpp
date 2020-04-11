#ifndef SRC_GHOST_OPERATOR_HPP_
#define SRC_GHOST_OPERATOR_HPP_

#include "grid.hpp"
#include "murphy.hpp"
#include "gridblock.hpp"

template <typename T>
using op_m_t = void (*)(const qid_t* qid, GridBlock* block, T);

/**
 * @brief The actual implementation of the Ghost iteration
 * 
 * @tparam O the type of operator that is used, see the definition of op_t
 * @tparam T the type of data given as user-defined data.
 * @tparam F the parameter pack containing the list of the fields which is taken as input/output, see the definition of opt_t
 * @param op the operator to be called on the ghosts
 * @param grid the grid which contains the ghosts
 * @param field the field on which we have to work
 * @param data the user-defined data forwared to the operator
 */
template <typename O, typename T>
void DoOp_M_(const O op, Grid* grid, T data) {
    m_begin;
    m_assert(grid->is_mesh_valid(), "mesh is not valid, unable to process");
    //-------------------------------------------------------------------------
    // get the grid info
    p8est_t*       forest  = grid->forest();
    p8est_ghost_t* ghost   = grid->ghost();
    const lid_t    nqlocal = ghost->mirrors.elem_count;  //number of ghost blocks

#pragma omp parallel for
    for (lid_t bid = 0; bid < nqlocal; bid++) {
        // get the mirror quad
        p8est_quadrant_t* mirror = p8est_quadrant_array_index(&ghost->mirrors, bid);
        // get the block and the tree
        p8est_tree_t*     tree   = p8est_tree_array_index(forest->trees, mirror->p.piggy3.which_tree);
        GridBlock*            block = reinterpret_cast<GridBlock*>(quad->p.user_data);
        // get the id, in this case the cummulative id = mirror id
        qid_t myid;
        myid.cid = bid;
        myid.qid = mirror->p.piggy3.local_num - tree->quadrants_offset;
        myid.tid = mirror->p.piggy3.which_tree;
        // send the task
        op(&myid, block, data);
    }
    //-------------------------------------------------------------------------
    m_end;
}

#endif  // SRC_GHOST_OPERATOR_HPP_