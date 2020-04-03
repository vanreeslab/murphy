#ifndef SRC_ITERATOR_HPP_
#define SRC_ITERATOR_HPP_

#include "block.hpp"
#include "field.hpp"
#include "grid.hpp"
#include "murphy.hpp"

template <class T>
using op_t = void (*)(Block* block, const qid_t* qid, const Field* fid, T);
template <class T>
using const_op_t = void (*)(const Block* block, const qid_t* qid, const Field* fid, T);
template <class T>
using bop_t = void (Block::*)(const qid_t* qid, Field* fid, T);

/**
 * @brief iterates over all the blocks in the forest and call the operator op
 * 
 * Because the operator is asusmed to change the field (and hence the block), we reset the ghost status
 * 
 * @tparam T the type of data which is given as context to the operator
 * @param op the operator itself, see definition of op_t
 * @param grid the grid on which one want to iterate
 * @param field the field on which the operation has to be done
 * @param data the context, i.e. user custom data to forward to the operator
 */
template <class T>
void Iterator(const op_t<T> op, Grid* grid, Field* field, T data) {
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
        Block*            block = (Block*)quad->p.user_data;
        // send the task
        op(block, &myid, field, data);
    }
    // downgrade the ghost status since we changed its value
    field->SetGhostStatus(false);
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief iterates over all the blocks in the forest and call the constant operator op
 * 
 * Because the operator is asusmed constant, i.e. it does NOT change the field, we don't change the ghost status
 * 
 * @tparam T the type of data which is given as context to the constant operator
 * @param op the operator itself, see definition of const_op_t
 * @param grid the grid on which one want to iterate
 * @param field the field on which the operation has to be done
 * @param data the context, i.e. user custom data to forward to the operator
 */
template <class T>
void Iterator(const const_op_t<T> op, Grid* grid, const Field* field, T data) {
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
        Block*            block = (Block*)quad->p.user_data;
        // send the task
        op(block, &myid, field, data);
    }
    //-------------------------------------------------------------------------
    m_end;
}

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
void Iterator(const bop_t<T> op, Grid* grid, Field* field, T data) {
    m_begin;
    //-------------------------------------------------------------------------
    // get the grid info
    p8est_t* forest = grid->forest();

    for (p4est_topidx_t it = forest->first_local_tree; it <= forest->last_local_tree; it++) {
        p8est_tree_t* tree    = p8est_tree_array_index(forest->trees, it);
        const lid_t   nqlocal = tree->quadrants.elem_count;
// #pragma omp parallel for
        for (lid_t bid = 0; bid < nqlocal; bid++) {
            // get the id
            qid_t myid;
            myid.cid = bid + tree->quadrants_offset;
            myid.qid = bid;
            myid.tid = it;
            // the quadrants can be from differents trees -> get the correct one
            p8est_quadrant_t* quad  = p8est_quadrant_array_index(&tree->quadrants, myid.qid);
            Block*            block = (Block*)quad->p.user_data;
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

#endif  // SRC_ITERATOR_HPP_
