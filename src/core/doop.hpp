#ifndef SRC_DOOP_HPP_
#define SRC_DOOP_HPP_

#include <type_traits>

#include "core/macros.hpp"
#include "core/types.hpp"
#include "grid/forestgrid.hpp"
#include "grid/gridblock.hpp"
#include "operator/blockoperator.hpp"
#include "toolsp4est.hpp"

/**
 * @brief iterates on the Gridblocks at a given level using the p4est_mesh object
 * 
 * This function is alike @ref DoOpMesh() but only calls the function on the given level
 * This functions takes an object @ref op and applies the function op->memfunc to it.
 * 
 * @tparam O the object type: std::nullptr_t if the object is nullptr or any pointer to an object
 * @tparam F the function type: a pointer to a member function of the GridBlock class (if @ref op is nullptr) or a pointer to a member function of the O class.
 * @tparam T the template package used to pass arguments to the functions
 * @param op a pointer to the object to use for the function call; if nullptr: the function has to be a member function of the GridBlock class
 * @param memfunc the member function to call with the GridBlock
 * @param grid the grid that contains the blocks
 * @param lvl the desired level
 * @param data the data passed by the user to the function
 */
template <typename O, typename F, typename... T>
void DoOpMesh(const O op, F memfunc, const ForestGrid*  grid, T... data) {
    m_begin;
    m_assert(grid->is_mesh_valid(), "mesh is not valid, unable to process");
    //-------------------------------------------------------------------------
    // do some static (=compilation) checks to be sure that the couple O and F is compatible
    constexpr bool do_gridblock = std::is_same<O, std::nullptr_t>();
    // if constexpr (do_gridblock) {
    //     constexpr bool is_member                = std::is_same<F, void (GridBlock::*)(T...)>();
    //     constexpr bool is_member_const          = std::is_same<F, void (GridBlock::*)(T...) const>();
    //     constexpr bool is_member_with_qid       = std::is_same<F, void (GridBlock::*)(const qid_t*, T...)>();
    //     constexpr bool is_member_const_with_qid = std::is_same<F, void (GridBlock::*)(const qid_t*, T...) const>();
    //     static_assert((is_member || is_member_const || is_member_with_qid || is_member_const_with_qid), "if the operator is nullptr, the function MUST be a member function of the GridBlock class");
    // } else {
    //     static_assert(std::is_pointer<O>(), "the operator type must be a pointer");
    //     static_assert(std::is_member_function_pointer<F>(), "the function type must be a pointer to a member function");
    //     constexpr bool is_member       = std::is_same<F, void (std::remove_pointer<O>::type::*)(const qid_t*, GridBlock*, T...)>();
    //     constexpr bool is_member_const = std::is_same<F, void (std::remove_pointer<O>::type::*)(const qid_t*, GridBlock*, T...) const>();
    //     static_assert(is_member || is_member_const, "if the operator is null, the function MUST be a member function of the GridBlock class");
    // }
    // check if we need to send the qid with it
    constexpr bool is_member_with_qid       = std::is_same<F, void (GridBlock::*)(const qid_t* , T...)>();
    constexpr bool is_member_const_with_qid = std::is_same<F, void (GridBlock::*)(const qid_t* , T...) const>();
    constexpr bool with_qid                 = is_member_const_with_qid || is_member_with_qid;
    //-------------------------------------------------------------------------
    // get the grid info
    p8est_t*      forest  = grid->p4est_forest();
    p8est_mesh_t* mesh    = grid->p4est_mesh();
    const iblock_t nqlocal = mesh->local_num_quadrants;  // number of trees * number of elem/tree

    for (iblock_t bid = 0; bid < nqlocal; ++bid) {
        // get the tree
        p8est_tree_t* tree = p8est_tree_array_index(forest->trees, mesh->quad_to_tree[bid]);
        // get the id
        qid_t myid;
        myid.cid = bid;                           // cummulative id
        myid.qid = bid - tree->quadrants_offset;  // quadrant id
        myid.tid = mesh->quad_to_tree[bid];       // tree id
        // the quadrants can be from differents trees -> get the correct one
        p8est_quadrant_t* quad;
        quad = p8est_quadrant_array_index(&tree->quadrants, myid.qid);

        GridBlock* block = *(reinterpret_cast<GridBlock**>(quad->p.user_data));
        // send the task on the block or on the operator, constexpr will compile only 1 of the two expressions
        if constexpr (do_gridblock && with_qid) {
            (block->*memfunc)(&myid, data...);
        } else if constexpr (do_gridblock) {
            (block->*memfunc)(data...);
        } else {
            (op->*memfunc)(&myid, block, data...);
        }
    }
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief iterates on the Gridblocks at a given level using the p4est_mesh object
 * 
 * This function is alike @ref DoOpMesh() but only calls the function on the given level
 * This functions takes an object @ref op and applies the function op->memfunc to it.
 * 
 * @tparam O the object type: std::nullptr_t if the object is nullptr or any pointer to an object
 * @tparam F the function type: a pointer to a member function of the GridBlock class (if @ref op is nullptr) or a pointer to a member function of the O class.
 * @tparam T the template package used to pass arguments to the functions
 * @param op a pointer to the object to use for the function call; if nullptr: the function has to be a member function of the GridBlock class
 * @param memfunc the member function to call with the GridBlock
 * @param grid the grid that contains the blocks
 * @param lvl the desired level
 * @param data the data passed by the user to the function
 */
template <typename O, typename F, typename... T>
void DoOpMeshLevel(const O op, F memfunc, const ForestGrid*  grid, const level_t lvl, T... data) {
    m_begin;
    m_assert(grid->is_mesh_valid(), "mesh is not valid, unable to process");
    //-------------------------------------------------------------------------
    // do some static (=compilation) checks to be sure that the couple O and F is compatible
    constexpr bool do_gridblock = std::is_same<O, std::nullptr_t>();
    // if constexpr (do_gridblock) {
    //     constexpr bool is_member                = std::is_same<F, void (GridBlock::*)(T...)>();
    //     constexpr bool is_member_const          = std::is_same<F, void (GridBlock::*)(T...) const>();
    //     constexpr bool is_member_with_qid       = std::is_same<F, void (GridBlock::*)(const qid_t*, T...)>();
    //     constexpr bool is_member_const_with_qid = std::is_same<F, void (GridBlock::*)(const qid_t*, T...) const>();
    //     static_assert((is_member || is_member_const || is_member_with_qid || is_member_const_with_qid), "if the operator is nullptr, the function MUST be a member function of the GridBlock class");
    // } else {
    //     static_assert(std::is_pointer<O>(), "the operator type must be a pointer");
    //     static_assert(std::is_member_function_pointer<F>(), "the function type must be a pointer to a member function");
    //     constexpr bool is_member       = std::is_same<F, void (std::remove_pointer<O>::type::*)(const qid_t*, GridBlock*, T...)>();
    //     constexpr bool is_member_const = std::is_same<F, void (std::remove_pointer<O>::type::*)(const qid_t*, GridBlock*, T...) const>();
    //     static_assert(is_member || is_member_const, "if the operator is null, the function MUST be a member function of the GridBlock class");
    // }
    // check if we need to send the qid with it
    constexpr bool is_member_with_qid       = std::is_same<F, void (GridBlock::*)(const qid_t* , T...)>();
    constexpr bool is_member_const_with_qid = std::is_same<F, void (GridBlock::*)(const qid_t* , T...) const>();
    constexpr bool with_qid                 = is_member_const_with_qid || is_member_with_qid;

    //-------------------------------------------------------------------------
    // get the grid info
    p8est_t*       forest  = grid->p4est_forest();
    p8est_mesh_t*  mesh    = grid->p4est_mesh();
    const iblock_t nqlocal = p4est_NumQuadOnLevel(mesh, lvl);

    for (iblock_t lid = 0; lid < nqlocal; ++lid) {
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

        m_assert(quad->level == lvl, "the selected quadrant has the wrong level");

        GridBlock* block = *(reinterpret_cast<GridBlock**>(quad->p.user_data));
        // send the task on the block or on the operator, constexpr will compile only 1 of the two expressions
        if constexpr (do_gridblock && with_qid) {
            (block->*memfunc)(&myid, data...);
        } else if constexpr (do_gridblock) {
            (block->*memfunc)(data...);
        } else {
            (op->*memfunc)(&myid, block, data...);
        }
    }
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief iterates on the Gridblocks using the p4est_tree, NOT using the p4est_mesh object
 * 
 * This function is alike @ref DoOpMesh() but can be called without the Ghost object created
 * 
 * @tparam O the object type: std::nullptr_t if the object is nullptr or any pointer to an object
 * @tparam F the function type: a pointer to a member function of the GridBlock class (if @ref op is nullptr) or a pointer to a member function of the O class.
 * @tparam T the template package used to pass arguments to the functions
 * @param op a pointer to the object to use for the function call; if nullptr: the function has to be a member function of the GridBlock class
 * @param memfunc the member function to call with the GridBlock
 * @param grid the grid that contains the blocks
 * @param data the data passed by the user to the function
 */
template <typename O, typename F, typename... T>
void DoOpTree(const O op, F memfunc, const ForestGrid*  grid, T... data) {
    m_begin;
    //-------------------------------------------------------------------------
    // do some static (=compilation) checks to be sure that the couple O and F is compatible
    constexpr bool do_gridblock = std::is_same<O, std::nullptr_t>();
    // if constexpr (do_gridblock) {
    //     constexpr bool is_member                = std::is_same<F, void (GridBlock::*)(T...)>();
    //     constexpr bool is_member_const          = std::is_same<F, void (GridBlock::*)(T...) const>();
    //     constexpr bool is_member_with_qid       = std::is_same<F, void (GridBlock::*)(const qid_t*, T...)>();
    //     constexpr bool is_member_const_with_qid = std::is_same<F, void (GridBlock::*)(const qid_t*, T...) const>();
    //     static_assert((is_member || is_member_const || is_member_with_qid || is_member_const_with_qid), "if the operator is nullptr, the function MUST be a member function of the GridBlock class");
    // } else {
    //     static_assert(std::is_pointer<O>(), "the operator type must be a pointer");
    //     static_assert(std::is_member_function_pointer<F>(), "the function type must be a pointer to a member function");
    //     constexpr bool is_member       = std::is_same<F, void (std::remove_pointer<O>::type::*)(const qid_t*, GridBlock*, T...)>();
    //     constexpr bool is_member_const = std::is_same<F, void (std::remove_pointer<O>::type::*)(const qid_t*, GridBlock*, T...) const>();
    //     static_assert(is_member || is_member_const, "if the operator is null, the function MUST be a member function of the GridBlock class");
    // }
    // check if we need to send the qid with it
    constexpr bool is_member_with_qid       = std::is_same<F, void (GridBlock::*)(const qid_t* , T...)>();
    constexpr bool is_member_const_with_qid = std::is_same<F, void (GridBlock::*)(const qid_t* , T...) const>();
    constexpr bool with_qid                 = is_member_const_with_qid || is_member_with_qid;
    //-------------------------------------------------------------------------
    // get the grid info
    p8est_t* forest = grid->p4est_forest();

    for (p4est_topidx_t it = forest->first_local_tree; it <= forest->last_local_tree; ++it) {
        p8est_tree_t* tree    = p8est_tree_array_index(forest->trees, it);
        const bidx_t  nqlocal = static_cast<bidx_t>(tree->quadrants.elem_count);

        for (bidx_t bid = 0; bid < nqlocal; ++bid) {
            p8est_quadrant_t* quad  = p8est_quadrant_array_index(&tree->quadrants, bid);
            GridBlock*        block = p4est_GetGridBlock(quad);

            // get the id
            qid_t myid;
            myid.cid = bid + tree->quadrants_offset;  // cummulative id
            myid.qid = bid;                           // quadrant id
            myid.tid = it;                            // tree id

            // send the task on the block or on the operator, constexpr will compile only 1 of the two expressions
            if constexpr (do_gridblock && with_qid) {
                (block->*memfunc)(&myid, data...);
            } else if constexpr (do_gridblock) {
                (block->*memfunc)(data...);
            } else {
                (op->*memfunc)(&myid, block, data...);
            }
        }
    }
    //-------------------------------------------------------------------------
    m_end;
}

#endif  // SRC_DOOP_HPP_