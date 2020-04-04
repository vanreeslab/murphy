
#ifndef SRC_OPERATOR_HPP_
#define SRC_OPERATOR_HPP_

#include "block.hpp"
#include "field.hpp"
#include "grid.hpp"
#include "murphy.hpp"

/**
 * @brief abstact class defining an operator
 * 
 */
class Operator {
   public:
    virtual void DoOp(const qid_t* qid, Block* block, Field* fid) = 0;
};

/**
 * @brief function used to perform the action of the Operator object.
 * 
 * For example, one calls an Operator using
 * ```
 * Operator* my_operator = new Operator();
 * DoOp<Operator>(CallOp,grid,field,my_operator);
 * ...
 * ```
 */
static void CallOp(const qid_t* qid, Block* block, Field* fid, Operator* op) {
    op->DoOp(qid, block, fid);
}

/**
 * @brief pointer to an operator function acting on a block
 * 
 * @tparam T the type of the user defined pointer used
 */
template <typename T>
using op_t = void (*)(const qid_t* qid, Block* block, Field* fid, T);
/**
 * @brief pointer to aconstant operator function acting on a block
 * 
 * @tparam T the type of the user defined pointer used
 */
template <typename T>
using const_op_t = void (*)(const qid_t* qid, const Block* block, const Field* fid, T);

/**
 * @brief implement the iteration on the Blocks
 * 
 * The iteration is performed using the mesh from the grid.
 * We use OpenMP to multithread the operation
 * 
 * @tparam O the type of operator considered (opt_t or const_op_t)
 * @tparam T the type of user-defined data fowarded to the operator function
 * @param op the actual operator function
 * @param grid the grid on which one want to iterate
 * @param field the field on which the operation has to be done
 * @param data the context, i.e. user custom data to forward to the operator
 */
template <typename O, typename T>
void DoOp_(const O op, Grid* grid, Field* field, T data) {
    m_begin;
    m_assert(grid->is_mesh_valid(), "mesh is not valid, unable to process");
    //-------------------------------------------------------------------------
    // get the grid info
    p8est_t*      forest  = grid->forest();
    p8est_mesh_t* mesh    = grid->mesh();
    const lid_t   nqlocal = mesh->local_num_quadrants;  // number of trees * number of elem/tree

#pragma omp parallel for
    for (lid_t bid = 0; bid < nqlocal; bid++) {
        // get the tree
        p8est_tree_t* tree = p8est_tree_array_index(forest->trees, mesh->quad_to_tree[bid]);
        // get the id
        qid_t myid;
        myid.cid = bid;
        myid.qid = bid - tree->quadrants_offset;
        myid.tid = mesh->quad_to_tree[bid];
        // the quadrants can be from differents trees -> get the correct one
        p8est_quadrant_t* quad  = p8est_quadrant_array_index(&tree->quadrants, myid.qid);
        Block*            block = reinterpret_cast<Block*>(quad->p.user_data);
        // send the task
        op(&myid, block, field, data);
    }
    // downgrade the ghost status if we are not a constant operator
    if (!std::is_same<O, const_op_t<T> >::value) {
        m_verb("setting the ghosts of %s to false", field->name().c_str());
        field->SetGhostStatus(false);
    }
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief wraps the function DoOp_() for operators, which performs the actual implementation
 * 
 * @tparam T the type of data which is given as context to the operator
 * @param op the operator itself, see definition of op_t
 * @param grid the grid on which one want to iterate
 * @param field the field on which the operation has to be done
 * @param data the context, i.e. user custom data to forward to the operator
 */
template <typename T>
void DoOp(const op_t<T> op, Grid* grid, Field* field, T data) {
    DoOp_<op_t<T>, T>(op, grid, field, data);
}

/**
 * @brief wraps the function DoOp_(), which performs the actual implementation
 * 
 * @tparam T the type of data which is given as context to the operator
 * @param op the constant operator itself, see definition of const_op_t
 * @param grid the grid on which one want to iterate
 * @param field the field on which the operation has to be done
 * @param data the context, i.e. user custom data to forward to the operator
 */
template <typename T>
void DoOp(const const_op_t<T> op, Grid* grid, Field* field, T data) {
    DoOp_<const_op_t<T>, T>(op, grid, field, data);
}

/**
 * @brief pointer to an member operator function of the class Block
 * 
 * @tparam T the type of the user-defined variable
 */
template <typename T>
using bop_t = void (Block::*)(const qid_t* qid, Field* fid, T);

/**
 * @brief iterates over all the blocks in the forest and call the block member fonction op
 * 
 * Because we operate on Block member functions, we assume that the mesh is not available.
 * Hence, we use the forest instead of the mesh to iterate on the blocks.
 * 
 * @tparam T the type of data which is given as context to the block function
 * @param op the operator itself, see definition of bop_t
 * @param grid the grid on which one want to iterate
 * @param field the field on which the operation has to be done
 * @param data the context, i.e. user custom data to forward to the operator
 */
template <class T>
void DoOp(const bop_t<T> op, Grid* grid, Field* field, T data) {
    m_begin;
    //-------------------------------------------------------------------------
    // get the grid info
    p8est_t* forest = grid->forest();

    for (p4est_topidx_t it = forest->first_local_tree; it <= forest->last_local_tree; it++) {
        p8est_tree_t* tree    = p8est_tree_array_index(forest->trees, it);
        const lid_t   nqlocal = tree->quadrants.elem_count;
#pragma omp parallel for
        for (lid_t bid = 0; bid < nqlocal; bid++) {
            // get the id
            qid_t myid;
            myid.cid = bid + tree->quadrants_offset;
            myid.qid = bid;
            myid.tid = it;
            // the quadrants can be from differents trees -> get the correct one
            p8est_quadrant_t* quad  = p8est_quadrant_array_index(&tree->quadrants, myid.qid);
            Block*            block = reinterpret_cast<Block*>(quad->p.user_data);
            // send the task
            // following https://en.cppreference.com/w/cpp/language/pointer
            (block->*op)(&myid, field, data);
        }
        // downgrade the ghost status since we changed its value
        field->SetGhostStatus(false);
        //-------------------------------------------------------------------------
        m_end;
    }
}

#endif  // SRC_OPERATOR_HPP_
