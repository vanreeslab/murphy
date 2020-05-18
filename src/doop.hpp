#ifndef SRC_DOOP_HPP_
#define SRC_DOOP_HPP_

#include "forestgrid.hpp"
#include "murphy.hpp"
#include "toolsp4est.hpp"
#include "gridblock.hpp"

/**
 * @brief General expression of any operator taking a list of fields as input/output
 * 
 * @warning This templated alias uses the parameter pack feature introduced in C++11
 * 
 * Here is some examples
 * - An operator `op_t<real_t*, Field*>` corresponds to 
 * ```
 * void operator1(const qid_t* qid, Block* block, Field* field1, real_t* data)
 * ```
 * 
 * - An operator `op_t<real_t*, const Field*, const Field*, Field*, Field*>` corresponds to
 * ```
 * void operator2(const qid_t* qid, Block* block, const Field* field_src1, const Field* field_src2, Field* field_trg1, Field* field_trg2, real_t* data)
 * ```
 * 
 * @tparam T the type of the data is passed as a user_defined pointer
 * @tparam F The parameter pack of Fields* which is passed
 */
template <typename T, typename... F>
using op_t = void (*)(const qid_t* qid, GridBlock* block, F... fid, T);

/**
 * @brief The actual implementation of the Block iteration
 * 
 * @warning We use a omp directive to loop over the blocks, hence the data is set as firstprivate!
 * 
 * @tparam O the type of operator that is used, see the definition of op_t
 * @tparam T the type of data given as user-defined data.
 * @tparam F the parameter pack containing the list of the fields which is taken as input/output, see the definition of opt_t
 * @param op the operator to be called on the blocks
 * @param grid the grid which contains the blocks
 * @param field the field on which we have to work
 * @param data the user-defined data forwared to the operator
 */
template <typename O, typename T, typename... F>
void DoOp_F_(const O op, ForestGrid* grid, F... field, T data) {
    m_begin;
    m_assert(grid->is_mesh_valid(), "mesh is not valid, unable to process");
    //-------------------------------------------------------------------------
    // get the grid info
    p8est_t*      forest  = grid->forest();
    p8est_mesh_t* mesh    = grid->mesh();
    const lid_t   nqlocal = mesh->local_num_quadrants;  // number of trees * number of elem/tree

#pragma omp parallel for firstprivate(data)
    for (lid_t bid = 0; bid < nqlocal; bid++) {
        // get the tree
        p8est_tree_t* tree;
        tree = p8est_tree_array_index(forest->trees, mesh->quad_to_tree[bid]);
        // get the id
        qid_t myid;
        myid.cid = bid;
        myid.qid = bid - tree->quadrants_offset;
        myid.tid = mesh->quad_to_tree[bid];
        // the quadrants can be from differents trees -> get the correct one
        p8est_quadrant_t* quad;
        quad = p8est_quadrant_array_index(&tree->quadrants, myid.qid);

        GridBlock* block = *(reinterpret_cast<GridBlock**>(quad->p.user_data));
        // send the task
        op(&myid, block, field..., data);
    }
    //-------------------------------------------------------------------------
    m_end;
}
/**
 * @brief The actual implementation of the Block iteration, performed ONLY on the given level
 * 
 * @warning We use a omp directive to loop over the blocks, hence the data is set as firstprivate!
 * 
 * @tparam O the type of operator that is used, see the definition of op_t
 * @tparam T the type of data given as user-defined data.
 * @tparam F the parameter pack containing the list of the fields which is taken as input/output, see the definition of opt_t
 * @param op the operator to be called on the blocks
 * @param grid the grid which contains the blocks
 * @param lvl the level on which we need to call the operator
 * @param field the field on which we have to work
 * @param data the user-defined data forwared to the operator
 */
template <typename O, typename T, typename... F>
void DoOp_F_(const O op, ForestGrid* grid, const level_t lvl, F... field, T data) {
    m_begin;
    m_assert(grid->is_mesh_valid(), "mesh is not valid, unable to process");
    //-------------------------------------------------------------------------
    // get the grid info
    p8est_t*       forest  = grid->forest();
    p8est_mesh_t*  mesh    = grid->mesh();
    const iblock_t nqlocal = p4est_NumQuadOnLevel(mesh, lvl);

#pragma omp parallel for firstprivate(data)
    for (iblock_t lid = 0; lid < nqlocal; lid++) {
        // get the corresponding id of the quadrant
        const iblock_t bid = p4est_GetQuadIdOnLevel(mesh, lvl, lid);
        // get the tree
        p8est_tree_t* tree;
        tree = p8est_tree_array_index(forest->trees, mesh->quad_to_tree[bid]);
        // get the id
        qid_t myid;
        myid.cid = bid;
        myid.qid = bid - tree->quadrants_offset;
        myid.tid = mesh->quad_to_tree[bid];
        // the quadrants can be from differents trees -> get the correct one
        p8est_quadrant_t* quad;
        quad = p8est_quadrant_array_index(&tree->quadrants, myid.qid);

        GridBlock* block = *(reinterpret_cast<GridBlock**>(quad->p.user_data));
        // send the task
        op(&myid, block, field..., data);
    }
    //-------------------------------------------------------------------------
    m_end;
}

#endif  // SRC_DOOP_HPP_