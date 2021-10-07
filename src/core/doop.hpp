#ifndef SRC_DOOP_HPP_
#define SRC_DOOP_HPP_

#include <type_traits>

#include "core/macros.hpp"
#include "core/types.hpp"
#include "grid/forestgrid.hpp"
#include "operator/blockoperator.hpp"
#include "tools/toolsp4est.hpp"
#include "tools/toolsblocktypes.hpp"

/**
 * @brief statically extract the expected block type from a member function pointer
 * 
 */
template<typename F> struct MemberBlockType;

template<typename O, typename BT, typename... Args>
struct MemberBlockType<void(O::*)(const qid_t*, BT*, Args...)>
{
    using BlockType = BT; 
};

template<typename O, typename BT, typename... Args>
struct MemberBlockType<void(O::*)(const qid_t*, BT*, Args...) const>
{
    using BlockType = BT; 
};

/**
 * @brief casts p4est user data to the correct pointer type before passing it to a member function
 * 
 * Does not perform any type checking, see CheckBlockType.
 * 
 * @tparam O class containing the member function
 * @tparam F type of the member function
 * @tparam T the template package used to pass arguments to the member function
 * @param op a pointer to the object to use for the function call
 * @param memfunc the member function to call on the block 
 * @param myid pointer to quadrant id for a given block
 * @param user_data void pointer to a given block
 * @param args the data passed by the user to the function
 */
template<typename O, typename F, typename... T>
void CallMemfunc(O op, F memfunc, qid_t* myid, void* user_data, T... args) {
    using BlockType = typename MemberBlockType<F>::BlockType;
    BlockType** p4est_usr_data = static_cast<BlockType**>(user_data);
    BlockType*  block          = p4est_usr_data[0];
    (op->*memfunc)(myid, block, args...);
}

/**
 * @brief checks at runtime that a member function is callable on the blocks in a Grid. 
 * 
 * @tparam F type of the member function
 * @param grid a pointer to the Grid
 * @param memfunc the member function to call on the block 
 */
template<typename F>
void CheckBlockType(const ForestGrid* grid, F memfunc) {
    using BlockType = typename MemberBlockType<F>::BlockType;
    const BlockDataType grid_block_type = grid->block_type();
    const BlockDataType func_block_type = TypeToEnum<BlockType>();
    m_assert(IsCompatibleBlockType(func_block_type, grid_block_type), 
        "argument and grid block types must be compatible");
}

//------------------------------------------------------------------------------
// defines the three different types of functions that can be called with DoOp

/**
 * @brief defines a function which is a member of the BlockType B
 */
template <typename B, typename... T>
using DoOpMesh_FnBlockMember = void (B::*)(T...);
/**
 * @brief defines a function which is a member of the BlockType B and requires the block id
 */
template <typename B, typename... T>
using DoOpMesh_FnBlockMemberQid = void (B::*)(const qid_t*, T...);
/**
 * @brief defines a function that needs a block id and the block B
 */
template <typename B, typename... T>
using DoOpMesh_FnOnBlockQid = void (*)(const qid_t*, B*, T...);

//------------------------------------------------------------------------------
/**
 * @brief call the function from the operator O on the block attachated to the quad
 */
template <typename O, typename F, typename... T>
void CallOperator(const BlockDataType grid_block_type,
                  const O op, F func, const qid_t* qid, const p8est_quadrant_t* quad, T... data) {
    m_assert(false, "This function is not callable through the DoOpMesh familly");
}

/**
 * @brief implements CallOperator for a DoOpMesh_FnBlockMember function
 */
template <typename B, typename... T>
void CallOperator(const BlockDataType  grid_block_type,
                  const std::nullptr_t op, DoOpMesh_FnBlockMember<B, T...>* func, const qid_t* qid, const p8est_quadrant_t* quad, T... data) {
    m_begin;
    m_assert(op == nullptr, "the operator must be nullptr");
    m_assert(IsCompatibleBlockType(TypeToEnum<B>(), grid_block_type), "The two types must be compatible to cast");
    //--------------------------------------------------------------------------
    B* block = p4est_GetBlock<B>(quad);
    block->func(data...);
    //--------------------------------------------------------------------------
    m_end;
}

/**
 * @brief implements CallOperator for a DoOpMesh_FnBlockMemberQid function
 */
template <typename B, typename... T>
void CallOperator(const BlockDataType  grid_block_type,
                  const std::nullptr_t op, DoOpMesh_FnBlockMemberQid<B, T...>* func, const qid_t* qid, const p8est_quadrant_t* quad, T... data) {
    m_begin;
    m_assert(op == nullptr, "the operator must be nullptr");
    m_assert(IsCompatibleBlockType(TypeToEnum<B>(), grid_block_type), "The two types must be compatible to cast");
    //--------------------------------------------------------------------------
    B* block = p4est_GetBlock<B>(quad);
    block->func(qid, data...);
    //--------------------------------------------------------------------------
    m_end;
}

/**
 * @brief implements CallOperator for a DoOpMesh_FnOnBlockQid function
 */
template <typename O, typename B, typename... T>
void CallOperator(const BlockDataType grid_block_type,
                  const O op, DoOpMesh_FnOnBlockQid<B, T...>* func, const qid_t* qid, const p8est_quadrant_t* quad, T... data) {
    m_begin;
    m_assert(op != nullptr, "the operator must be nullptr");
    m_assert(IsCompatibleBlockType(TypeToEnum<B>(), grid_block_type), "The two types must be compatible to cast");
    //--------------------------------------------------------------------------
    B* block = p4est_GetBlock<B>(quad);
    op->func(qid, block, data...);
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
    //--------------------------------------------------------------------------
    // // do some static (=compilation) checks to be sure that the couple O and F is compatible
    // constexpr bool do_gridblock = std::is_same<O, std::nullptr_t>();
    // // if constexpr (do_gridblock) {
    // //     constexpr bool is_member                = std::is_same<F, void (GridBlock::*)(T...)>();
    // //     constexpr bool is_member_const          = std::is_same<F, void (GridBlock::*)(T...) const>();
    // //     constexpr bool is_member_with_qid       = std::is_same<F, void (GridBlock::*)(const qid_t*, T...)>();
    // //     constexpr bool is_member_const_with_qid = std::is_same<F, void (GridBlock::*)(const qid_t*, T...) const>();
    // //     static_assert((is_member || is_member_const || is_member_with_qid || is_member_const_with_qid), "if the operator is nullptr, the function MUST be a member function of the GridBlock class");
    // // } else {
    // //     static_assert(std::is_pointer<O>(), "the operator type must be a pointer");
    // //     static_assert(std::is_member_function_pointer<F>(), "the function type must be a pointer to a member function");
    // //     constexpr bool is_member       = std::is_same<F, void (std::remove_pointer<O>::type::*)(const qid_t*, GridBlock*, T...)>();
    // //     constexpr bool is_member_const = std::is_same<F, void (std::remove_pointer<O>::type::*)(const qid_t*, GridBlock*, T...) const>();
    // //     static_assert(is_member || is_member_const, "if the operator is null, the function MUST be a member function of the GridBlock class");
    // // }
    // // check if we need to send the qid with it
    // constexpr bool is_member_with_qid       = std::is_same<F, void (GridBlock::*)(const qid_t* , T...)>();
    // constexpr bool is_member_const_with_qid = std::is_same<F, void (GridBlock::*)(const qid_t* , T...) const>();
    // constexpr bool with_qid                 = is_member_const_with_qid || is_member_with_qid;
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

        // // send the task on the block or on the operator, constexpr will compile only 1 of the two expressions
        // if constexpr (do_gridblock && with_qid) {
        //     GridBlock* block = p4est_GetBlock(quad);
        //     (block->*memfunc)(&myid, data...);
        // } else if constexpr (do_gridblock) {
        //     GridBlock* block = p4est_GetBlock(quad);
        //     (block->*memfunc)(data...);
        // } else {
        //     CheckBlockType(grid, memfunc);
        //     CallMemfunc(op, memfunc, &myid, quad->p.user_data, data...);
        // }
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
    // // do some static (=compilation) checks to be sure that the couple O and F is compatible
    // constexpr bool do_gridblock = std::is_same<O, std::nullptr_t>();
    // // if constexpr (do_gridblock) {
    // //     constexpr bool is_member                = std::is_same<F, void (GridBlock::*)(T...)>();
    // //     constexpr bool is_member_const          = std::is_same<F, void (GridBlock::*)(T...) const>();
    // //     constexpr bool is_member_with_qid       = std::is_same<F, void (GridBlock::*)(const qid_t*, T...)>();
    // //     constexpr bool is_member_const_with_qid = std::is_same<F, void (GridBlock::*)(const qid_t*, T...) const>();
    // //     static_assert((is_member || is_member_const || is_member_with_qid || is_member_const_with_qid), "if the operator is nullptr, the function MUST be a member function of the GridBlock class");
    // // } else {
    // //     static_assert(std::is_pointer<O>(), "the operator type must be a pointer");
    // //     static_assert(std::is_member_function_pointer<F>(), "the function type must be a pointer to a member function");
    // //     constexpr bool is_member       = std::is_same<F, void (std::remove_pointer<O>::type::*)(const qid_t*, GridBlock*, T...)>();
    // //     constexpr bool is_member_const = std::is_same<F, void (std::remove_pointer<O>::type::*)(const qid_t*, GridBlock*, T...) const>();
    // //     static_assert(is_member || is_member_const, "if the operator is null, the function MUST be a member function of the GridBlock class");
    // // }
    // // check if we need to send the qid with it
    // constexpr bool is_member_with_qid       = std::is_same<F, void (GridBlock::*)(const qid_t* , T...)>();
    // constexpr bool is_member_const_with_qid = std::is_same<F, void (GridBlock::*)(const qid_t* , T...) const>();
    // constexpr bool with_qid                 = is_member_const_with_qid || is_member_with_qid;

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

        // // send the task on the block or on the operator, constexpr will compile only 1 of the two expressions
        // if constexpr (do_gridblock && with_qid) {
        //     GridBlock* block = p4est_GetBlock(quad);
        //     (block->*memfunc)(&myid, data...);
        // } else if constexpr (do_gridblock) {
        //     GridBlock* block = p4est_GetBlock(quad);
        //     (block->*memfunc)(data...);
        // } else {
        //     CheckBlockType(grid, memfunc);
        //     CallMemfunc(op, memfunc, &myid, quad->p.user_data, data...);
        // }
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

            // get the id
            qid_t myid;
            myid.cid = bid + tree->quadrants_offset;  // cummulative id
            myid.qid = bid;                           // quadrant id
            myid.tid = it;                            // tree id

            // // send the task on the block or on the operator, constexpr will compile only 1 of the two expressions
            // if constexpr (do_gridblock && with_qid) {
            //     GridBlock* block = p4est_GetBlock<GridBlock>(quad);
            //     (block->*memfunc)(&myid, data...);
            // } else if constexpr (do_gridblock) {
            //     GridBlock* block = p4est_GetBlock<GridBlock>(quad);
            //     (block->*memfunc)(data...);
            // } else {
            //     CheckBlockType(grid, memfunc);
            //     CallMemfunc(op, memfunc, &myid, quad->p.user_data, data...);
            // }
            CallOperator(grid->block_type(), op, memfunc, &myid, quad, data...);
        }
    }
    //-------------------------------------------------------------------------
    m_end;
}

#endif  // SRC_DOOP_HPP_