#ifndef SRC_DOOP_HPP_
#define SRC_DOOP_HPP_

#include <type_traits>

#include "core/macros.hpp"
#include "core/types.hpp"
#include "grid/forestgrid.hpp"
#include "operator/blockoperator.hpp"
#include "tools/toolsp4est.hpp"
#include "tools/toolsblocktypes.hpp"

//------------------------------------------------------------------------------
// defines the three different types of functions that can be called with DoOp

/**
 * @brief defines a function which is a member of the BlockType B
 */
template <typename B, typename... T>
using DoOpMesh_FnBlockMember = void (B::*)(T...);
template <typename B, typename... T>
using DoOpMesh_FnBlockMemberConst = void (B::*)(T...) const;

/**
 * @brief defines a function which is a member of the BlockType B and requires the block id
 */
template <typename B, typename... T>
using DoOpMesh_FnBlockMemberQid = void (B::*)(const qid_t*, T...);
template <typename B, typename... T>
using DoOpMesh_FnBlockMemberQidConst = void (B::*)(const qid_t*, T...) const;

/**
 * @brief defines a function that needs a block id and the block B
 * 
 * @note std::remove_pointer<O> is required here!! Example: O is set to SetValue* and we want a function on SetValue
 */
template <typename O, typename B, typename... T>
using DoOpMesh_FnOnBlockQid = void (std::remove_pointer<O>::type::*)(const qid_t*, B*, T...);
template <typename O, typename B, typename... T>
using DoOpMesh_FnOnBlockQidConst = void (std::remove_pointer<O>::type::*)(const qid_t*, B*, T...)const;

//------------------------------------------------------------------------------
/**
 * @brief call the function from the operator O on the block attachated to the quad
 */
template <typename O, typename F, typename... TD>
void CallOperator(const BlockDataType grid_block_type,
                  const O op, F func, const qid_t* qid, const p8est_quadrant_t* quad, TD... data) {
    m_assert(false, "This function < %s > is not callable in %s",typeid(func).name(),__PRETTY_FUNCTION__);
}

/**
 * @brief implements CallOperator for a DoOpMesh_FnBlockMember function
 */
template <typename B, typename... TF, typename... TD>
void CallOperator(const BlockDataType  grid_block_type,
                  const std::nullptr_t op, DoOpMesh_FnBlockMember<B, TF...> func, const qid_t* qid, const p8est_quadrant_t* quad, TD... data) {
    m_begin;
    m_assert(op == nullptr, "the operator must be nullptr");
    m_assert(IsCompatibleBlockType(TypeToEnum<B>(), grid_block_type), "The two types must be compatible to cast");
    static_assert(sizeof...(TD) ==  sizeof...(TF), "The size of the arguments for the DoOpMesh function and for the operator must be the same");
    //--------------------------------------------------------------------------
    B* block = p4est_GetBlock<B>(quad);
    (block->*func)(data...);
    //--------------------------------------------------------------------------
    m_end;
}
/**
 * @brief implements CallOperator for a DoOpMesh_FnBlockMemberConst function
 */
template <typename B, typename... TF, typename... TD>
void CallOperator(const BlockDataType  grid_block_type,
                  const std::nullptr_t op, DoOpMesh_FnBlockMemberConst<B, TF...> func, const qid_t* qid, const p8est_quadrant_t* quad, TD... data) {
    m_begin;
    m_assert(op == nullptr, "the operator must be nullptr");
    m_assert(IsCompatibleBlockType(TypeToEnum<B>(), grid_block_type), "The two types must be compatible to cast");
    static_assert(sizeof...(TD) ==  sizeof...(TF), "The size of the arguments for the DoOpMesh function and for the operator must be the same");
    //--------------------------------------------------------------------------
    B* block = p4est_GetBlock<B>(quad);
    (block->*func)(data...);
    //--------------------------------------------------------------------------
    m_end;
}

/**
 * @brief implements CallOperator for a DoOpMesh_FnBlockMemberQid function
 */
template <typename B, typename... TF, typename... TD>
void CallOperator(const BlockDataType  grid_block_type,
                  const std::nullptr_t op, DoOpMesh_FnBlockMemberQid<B, TF...> func, const qid_t* qid, const p8est_quadrant_t* quad, TD... data) {
    m_begin;
    m_assert(op == nullptr, "the operator must be nullptr");
    m_assert(IsCompatibleBlockType(TypeToEnum<B>(), grid_block_type), "The two types must be compatible to cast");
    static_assert(sizeof...(TD) ==  sizeof...(TF), "The size of the arguments for the DoOpMesh function and for the operator must be the same");
    //--------------------------------------------------------------------------
    B* block = p4est_GetBlock<B>(quad);
    (block->*func)(qid, data...);
    //--------------------------------------------------------------------------
    m_end;
}
/**
 * @brief implements CallOperator for a DoOpMesh_FnBlockMemberQidConst function
 */
template <typename B, typename... TF, typename... TD>
void CallOperator(const BlockDataType  grid_block_type,
                  const std::nullptr_t op, DoOpMesh_FnBlockMemberQidConst<B, TF...> func, const qid_t* qid, const p8est_quadrant_t* quad, TD... data) {
    m_begin;
    m_assert(op == nullptr, "the operator must be nullptr");
    m_assert(IsCompatibleBlockType(TypeToEnum<B>(), grid_block_type), "The two types must be compatible to cast");
    static_assert(sizeof...(TD) ==  sizeof...(TF), "The size of the arguments for the DoOpMesh function and for the operator must be the same");
    //--------------------------------------------------------------------------
    B* block = p4est_GetBlock<B>(quad);
    (block->*func)(qid, data...);
    //--------------------------------------------------------------------------
    m_end;
}

/**
 * @brief implements CallOperator for a DoOpMesh_FnOnBlockQid function
 */
template <typename O, typename B, typename... TF, typename... TD>
void CallOperator(const BlockDataType grid_block_type,
                  const O op, DoOpMesh_FnOnBlockQid<O, B, TF...> func, const qid_t* qid, const p8est_quadrant_t* quad, TD... data) {
    m_begin;
    m_assert(op != nullptr, "the operator must be nullptr");
    m_assert(IsCompatibleBlockType(TypeToEnum<B>(), grid_block_type), "The two types must be compatible to cast");
    static_assert(sizeof...(TD) ==  sizeof...(TF), "The size of the arguments for the DoOpMesh function and for the operator must be the same");
    //--------------------------------------------------------------------------
    B* block = p4est_GetBlock<B>(quad);
    (op->*func)(qid, block, data...);
    //--------------------------------------------------------------------------
    m_end;
}
/**
 * @brief implements CallOperator for a DoOpMesh_FnOnBlockQid function
 */
template <typename O, typename B, typename... TF, typename... TD>
void CallOperator(const BlockDataType grid_block_type,
                  const O op, DoOpMesh_FnOnBlockQidConst<O, B, TF...> func, const qid_t* qid, const p8est_quadrant_t* quad, TD... data) {
    m_begin;
    m_assert(op != nullptr, "the operator must be nullptr");
    m_assert(IsCompatibleBlockType(TypeToEnum<B>(), grid_block_type), "The two types must be compatible to cast");
    static_assert(sizeof...(TD) ==  sizeof...(TF), "The size of the arguments for the DoOpMesh function and for the operator must be the same");
    //--------------------------------------------------------------------------
    B* block = p4est_GetBlock<B>(quad);
    (op->*func)(qid, block, data...);
    //--------------------------------------------------------------------------
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
void DoOpMesh(const O op, F memfunc, const ForestGrid*  grid, T... data) {
    m_begin;
    m_assert(grid->is_mesh_valid(), "mesh is not valid, unable to process");
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
        p8est_quadrant_t* quad = p8est_quadrant_array_index(&tree->quadrants, myid.qid);
        CallOperator(grid->block_type(), op, memfunc, &myid, quad, data...);
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
        p8est_quadrant_t* quad  = p8est_quadrant_array_index(&tree->quadrants, myid.qid);
        m_assert(quad->level == lvl, "the selected quadrant has the wrong level");

        // send the task on the block or on the operator
        CallOperator(grid->block_type(), op, memfunc, &myid, quad, data...);
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
    // get the grid info
    p8est_t* forest = grid->p4est_forest();

    for (p4est_topidx_t it = forest->first_local_tree; it <= forest->last_local_tree; ++it) {
        p8est_tree_t* tree    = p8est_tree_array_index(forest->trees, it);
        const bidx_t  nqlocal = static_cast<bidx_t>(tree->quadrants.elem_count);

        for (bidx_t bid = 0; bid < nqlocal; ++bid) {
            p8est_quadrant_t* quad  = p8est_quadrant_array_index(&tree->quadrants, bid);

            // get the id
            qid_t myid;
            myid.cid = bid + tree->quadrants_offset;  // cummulative id
            myid.qid = bid;                           // quadrant id
            myid.tid = it;                            // tree id
            
            CallOperator(grid->block_type(), op, memfunc, &myid, quad, data...);
        }
    }
    //-------------------------------------------------------------------------
    m_end;
}

#endif  // SRC_DOOP_HPP_