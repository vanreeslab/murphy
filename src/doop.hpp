#ifndef SRC_DOOP_HPP_
#define SRC_DOOP_HPP_

#include "forestgrid.hpp"
#include "murphy.hpp"
#include "toolsp4est.hpp"
#include "gridblock.hpp"

/**
 * @brief General expression of any operator taking a list of types as input
 * 
 * @warning This templated alias uses the parameter pack feature introduced in C++11
 * 
 * Here is some examples
 * - An operator `op_t<Field*, real_t*>` corresponds to 
 * ```
 * void operator1(const qid_t* qid, Block* block, Field* field1, real_t* data)
 * ```
 * 
 * - An operator `op_t<const Field*, const Field*, Field*, Field*, real_t*>` corresponds to
 * ```
 * void operator2(const qid_t* qid, Block* block, const Field* field_src1, const Field* field_src2, Field* field_trg1, Field* field_trg2, real_t* data)
 * ```
 * 
 * @tparam T the type of the data is passed after the qid and the GridBlock
 */
template <typename... T>
using op_t = void (*)(const qid_t* qid, GridBlock* block, T... data);

/**
 * @brief Iterates on every GridBlock using the p4est_mesh
 * 
 * @tparam T the parameter pack containing the list of the data to transfer to the operator
 * @param op the operator to be called on the blocks
 * @param grid the grid which contains the blocks
 * @param data the data forwared to the operator
 */
template <typename... T>
void DoOpMesh(const op_t<T...> op, ForestGrid* grid, T... data) {
    m_begin;
    m_assert(grid->is_mesh_valid(), "mesh is not valid, unable to process");
    //-------------------------------------------------------------------------
    // get the grid info
    p8est_t*      forest  = grid->p4est_forest();
    p8est_mesh_t* mesh    = grid->p4est_mesh();
    const lid_t   nqlocal = mesh->local_num_quadrants;  // number of trees * number of elem/tree

    for (lid_t bid = 0; bid < nqlocal; bid++) {
        // get the tree
        p8est_tree_t* tree;
        tree = p8est_tree_array_index(forest->trees, mesh->quad_to_tree[bid]);
        // get the id
        qid_t myid;
        myid.cid = bid;                           // cummulative id
        myid.qid = bid - tree->quadrants_offset;  // quadrant id
        myid.tid = mesh->quad_to_tree[bid];       // tree id
        // the quadrants can be from differents trees -> get the correct one
        p8est_quadrant_t* quad;
        quad = p8est_quadrant_array_index(&tree->quadrants, myid.qid);

        GridBlock* block = *(reinterpret_cast<GridBlock**>(quad->p.user_data));
        // send the task
        (*op)(&myid, block, data...);
    }
    //-------------------------------------------------------------------------
    m_end;
}
/**
 * @brief Iterates on GridBlock at a given level using the p4est_mesh
 * 
 * @warning We use a omp directive to loop over the blocks, hence the data is set as firstprivate!
 * 
 * @tparam T the parameter pack containing the list of the data to transfer to the operator
 * @param op the operator to be called on the blocks
 * @param grid the grid which contains the blocks
 * @param lvl the level on which we need to call the operator
 * @param data the data forwared to the operator
 */
template <typename... T>
void DoOpMeshLevel(const op_t<T...> op, ForestGrid* grid, const level_t lvl, T... data) {
    m_begin;
    m_assert(grid->is_mesh_valid(), "mesh is not valid, unable to process");
    //-------------------------------------------------------------------------
    // get the grid info
    p8est_t*       forest  = grid->p4est_forest();
    p8est_mesh_t*  mesh    = grid->p4est_mesh();
    const iblock_t nqlocal = p4est_NumQuadOnLevel(mesh, lvl);

    for (iblock_t lid = 0; lid < nqlocal; lid++) {
        // get the corresponding id of the quadrant
        const iblock_t bid = p4est_GetQuadIdOnLevel(mesh, lvl, lid);
        // get the tree
        p8est_tree_t* tree = p8est_tree_array_index(forest->trees, mesh->quad_to_tree[bid]);
        // get the id
        qid_t myid;
        myid.cid = bid;
        myid.qid = bid - tree->quadrants_offset;
        myid.tid = mesh->quad_to_tree[bid];
        // the quadrants can be from differents trees -> get the correct one
        p8est_quadrant_t* quad = p8est_quadrant_array_index(&tree->quadrants, myid.qid);

        GridBlock* block = *(reinterpret_cast<GridBlock**>(quad->p.user_data));
        // send the task
        (*op)(&myid, block, data...);
    }
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief Iterates on every GridBlock NOT using the p4est_mesh
 * 
 * @tparam T the parameter pack containing the list of the data to transfer to the operator
 * @param op the operator to be called on the blocks
 * @param grid the grid which contains the blocks
 * @param data the data forwared to the operator
 */
template <typename... T>
void DoOpTree(const op_t<T...> op, ForestGrid* grid, T... data) {
    m_begin;
    //-------------------------------------------------------------------------
    // get the grid info
    p8est_t* forest = grid->p4est_forest();

    for (p4est_topidx_t it = forest->first_local_tree; it <= forest->last_local_tree; ++it) {
        p8est_tree_t* tree    = p8est_tree_array_index(forest->trees, it);
        const size_t  nqlocal = tree->quadrants.elem_count;

        for (size_t bid = 0; bid < nqlocal; bid++) {
            p8est_quadrant_t* quad  = p8est_quadrant_array_index(&tree->quadrants, bid);
            GridBlock*        block = *(reinterpret_cast<GridBlock**>(quad->p.user_data));

            // get the id
            qid_t myid;
            myid.cid = bid;                           // cummulative id
            myid.qid = bid + tree->quadrants_offset;  // quadrant id
            myid.tid = it;                            // tree id
            // send the task
            (*op)(&myid, block, data...);
        }
        //-------------------------------------------------------------------------
        m_end;
    }
}

#endif  // SRC_DOOP_HPP_