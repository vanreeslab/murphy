
/**
 * @file operator.hpp
 * 
 * @warning This is advanced feature for Murphy that rely on parameter packs (C++11 feature).
 * Please consider any change to this file very carefully
 * 
 * There are two foreseen, yet equivalent, ways to interact with blocks using `Operators`:
 * - create an object that inherites from one of the available class (OperatorF, ConstOperatorF and OperatorF2F)
 * - create a function that satisfies the requirement imposed (op_t<T, Field*>, op_t<T, const Field*> and op_t<T,const Field*,Field*>)
 * 
 * In both approaches, a chain of template will be used to redirect the requirement to the function DoOp_().
 * 
 * We highly recommand the user to take the object approach.
 * In this way of thincking, the user inheritates form one of the three objects: OperatorF, ConstOperatorF and OperatorF2F.
 * One has just to provide an implementation of the function apply().
 * 
 * However, this is not always possible or convenient, hence the function-based approach has been preserved.
 * More information can be found in DoOp() functions and @ref op_t definition
 * 
 * Finally, a last way has been implemented, restricted for the Block member function.
 * We strongly advise the user to **NOT** use is as it relies on a less efficient looping strategies, specific for the use of the Block class.
 * 
 * 
 */

#ifndef SRC_OPERATOR_HPP_
#define SRC_OPERATOR_HPP_

#include "block.hpp"
#include "field.hpp"
#include "grid.hpp"
#include "murphy.hpp"

//=================================================================================================
/**
 * @brief A simple operator on a block, which modifies its content
 * 
 */
class OperatorF {
   public:
   /**
    * @brief Implementation of this virtual function has to be provided by the user as a member function
    */
    virtual void apply(const qid_t* qid, Block* block, Field* fid) = 0;
    /**
     * @brief call OperatorF::apply() on each block and change the ghost status of Field to `false`
     */
    void         DoOp(Grid* grid, Field* field);
};
/**
 * @brief this function is called by DoOp_() function (through OperatorF::DoOp()) to apply the operation to a considered Block
 * 
 * @param qid the reference of the block, see qid_t
 * @param block the Block itself
 * @param fid the field on which we operate
 * @param op the OperatorF object containing all the needed data
 */
void CallOpF(const qid_t* qid, Block* block, Field* fid, OperatorF* op);

//=================================================================================================
/**
 * @brief a constant Operator, i.e. which does not modify the Field
 * 
 */
class ConstOperatorF {
   public:
   /**
    * @brief Implementation of this virtual function has to be provided by the user as a member function
    */
    virtual void apply(const qid_t* qid, Block* block, const Field* fid) = 0;
    /**
     * @brief call ConstOperatorF::apply() on each block
     */
    void         DoOp(Grid* grid, Field* field);
};
/**
 * @brief this function is called by DoOp_() function (through ConstOperatorF::DoOp()) to apply the operation to a considered Block
 * 
 * @param qid the reference of the block, see qid_t
 * @param block the Block itself, which cannot be modified
 * @param fid the field on which we operate
 * @param op the ConstOperatorF object containing all the needed data
 */
void ConstCallOpF(const qid_t* qid, Block* block, const Field* fid, ConstOperatorF* op);

//=================================================================================================
/**
 * @brief a Field to Field operator, i.e. which uses the content of one Field to modify another Field
 * 
 */
class OperatorF2F {
   public:
   /**
    * @brief Implementation of this virtual function has to be provided by the user as a member function
    */
    virtual void apply(const qid_t* qid, Block* block, const Field* fid_src, Field* fid_trg) = 0;
    /**
     * @brief call OperatorF2F::apply() on each block and change the ghost status of the Field field_trg to `false`
     */
    void         DoOp(Grid* grid, const Field* field_src, Field* field_trg);
};
/**
 * @brief this function is called by DoOp_() function (through OperatorF2F::DoOp()) to apply the operation to a considered Block
 * 
 * @param qid the reference of the block, see qid_t
 * @param block the Block itself, which cannot be modified
 * @param fid_src the source field
 * @param fid_trg the target field
 * @param op the OperatorF2F object containing all the needed data
 */
void CallOpF2F(const qid_t* qid, Block* block, const Field* fid_src, Field* fid_trg, OperatorF2F* op);



//=================================================================================================
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
using op_t = void (*)(const qid_t* qid, Block* block, F... fid, T);

/**
 * @brief The actual implementation of the Block iteration
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
void DoOp_(const O op, Grid* grid, F... field, T data) {
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
        op(&myid, block, field..., data);
    }
    //-------------------------------------------------------------------------
    m_end;
}

//=================================================================================================
/**
 * @brief wraps the function DoOp_() for an operator (`op_t<T, Field*>`, see @ref op_t), which performs the actual implementation
 * 
 * We first call the operation on the blocks, then set the ghosts as changed
 * 
 * @tparam T the type of data which is given as context to the operator
 * @param op the operator itself, see definition of op_t
 * @param grid the grid on which one want to iterate
 * @param field the field on which the operation has to be done
 * @param data the context, i.e. user custom data to forward to the operator
 */
template <typename T>
void DoOp(const op_t<T, Field*> op, Grid* grid, Field* field, T data) {
    // execute the operator on the blocks
    DoOp_<op_t<T, Field*>, T, Field*>(op, grid, field, data);
    // set the ghost as changed
    m_verb("setting the ghosts of %s to false", field->name().c_str());
    field->SetGhostStatus(false);
}

/**
 * @brief wraps the function DoOp_() for a constant operator (`const op_t<T, const Field*>`, see @ref op_t), which performs the actual implementation
 * 
 * 
 * @tparam T the type of data which is given as context to the operator
 * @param op the operator itself, see definition of op_t
 * @param grid the grid on which one want to iterate
 * @param field the field on which the operation has to be done
 * @param data the context, i.e. user custom data to forward to the operator
 */
template <typename T>
void DoOp(const op_t<T, const Field*> op, Grid* grid, Field* field, T data) {
    // execute the operator on the blocks
    DoOp_<op_t<T, const Field*>, T, const Field*>(op, grid, field, data);
    // no change of the ghost is needed
}

/**
 * @brief wraps the function DoOp_() for a field to field operator (`op_t<T, const Field*, Field*>`, see @ref op_t), which performs the actual implementation
 * 
 * We first call the operation on the blocks, then set the ghosts as changed
 * 
 * @tparam T the type of data which is given as context to the operator
 * @param op the operator itself, see definition of op_t
 * @param grid the grid on which one want to iterate
 * @param field the field on which the operation has to be done
 * @param data the context, i.e. user custom data to forward to the operator
 */
template <typename T>
void DoOp(const op_t<T, const Field*, Field*> op, Grid* grid, Field* field_src, Field* field_trg, T data) {
    // execute the operation on the blocks
    DoOp_<op_t<T, const Field*, Field*>, T, const Field*, Field*>(op, grid, field_src, field_trg, data);
    // set the field_trg ghosts as changed
    m_verb("setting the ghosts of %s to false", field_trg->name().c_str());
    field_trg->SetGhostStatus(false);
}


//=========================================================================================
/**
 * @brief pointer to an member function of the class Block
 * 
 * @tparam T the type of the user-defined variable
 */
template <typename T>
using bop_t = void (Block::*)(const qid_t* qid, Field* fid, T);

/**
 * @brief iterates over all the blocks in the forest and call the block member function provided
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
            // send the task following https://en.cppreference.com/w/cpp/language/pointer
            (block->*op)(&myid, field, data);
        }
        // downgrade the ghost status since we changed its value
        field->SetGhostStatus(false);
        //-------------------------------------------------------------------------
        m_end;
    }
}

#endif  // SRC_OPERATOR_HPP_
