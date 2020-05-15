#ifndef SRC_OPERATOR_HPP_
#define SRC_OPERATOR_HPP_

#include <limits>

#include "field.hpp"
#include "forestgrid.hpp"
#include "gridblock.hpp"
#include "murphy.hpp"
#include "toolsp4est.hpp"

using std::nullptr_t;
using std::numeric_limits;

//=================================================================================================
/**
 * @brief A simple operator on a block, which does not use a field
 * 
 */
class OperatorS {
   public:
    /**
     * @brief Implementation of this virtual function has to be provided by the user as a member function
     * 
     * @warning this function is processed with a multi-thread environment
     * 
     * @param qid the id of the quadrant
     * @param block the block itself
     */
    virtual void ApplyOpS(const qid_t* qid, GridBlock* block) = 0;
    /**
     * @brief call OperatorF::ApplyOpS() on each block and change the ghost status of Field to `false`
     */
    virtual void operator()(ForestGrid* grid);
};
/**
 * @brief this function is called by DoOp_() function (through OperatorF::operator()()) to apply the operation to a considered Block
 * 
 * @param qid the reference of the block, see qid_t
 * @param block the Block itself
 * @param fid it's value has to be nullptr, still we keep it to match the @ref op_t definition
 * @param op the OperatorS object containing all the needed data
 */
void CallOpS(const qid_t* qid, GridBlock* block, nullptr_t* fid, OperatorS* op);

//=================================================================================================
/**
 * @brief A field operator on a block, which modifies its content
 * 
 */
class OperatorF {
   public:
    /**
     * @brief Implementation of this virtual function has to be provided by the user as a member function
     * 
     * @warning this function is processed with a multi-thread environment
     * 
     * @param qid the id of the quadrant
     * @param block the block itself
     * @param fid the field on which we execute the operation
     */
    virtual void ApplyOpF(const qid_t* qid, GridBlock* block, Field* fid) = 0;
    /**
     * @brief call OperatorF::ApplyOpF() on each block and change the ghost status of Field to `false`
     */
    virtual void operator()(ForestGrid* grid, Field* field);
};
/**
 * @brief this function is called by DoOp_() function (through OperatorF::operator()()) to apply the operation to a considered Block
 * 
 * @param qid the reference of the block, see qid_t
 * @param block the Block itself
 * @param fid the field on which we operate
 * @param op the OperatorF object containing all the needed data
 */
void CallOpF(const qid_t* qid, GridBlock* block, Field* fid, OperatorF* op);

//=================================================================================================
/**
 * @brief a constant field Operator, i.e. which does not modify the considered field
 * 
 */
class ConstOperatorF {
   public:
    /**
     * @brief Implementation of this virtual function has to be provided by the user as a member function
     * 
     * @warning this function is processed with a multi-thread environment
     * 
     * @param qid the id of the quadrant
     * @param block the block itself
     * @param fid the field on which we execute it
     */
    virtual void ApplyConstOpF(const qid_t* qid, GridBlock* block, const Field* fid) = 0;
    /**
     * @brief call ConstOperatorF::ApplyConstOpF() on each block
     */
    virtual void operator()(ForestGrid* grid, Field* field);
};
/**
 * @brief this function is called by DoOp_() function (through ConstOperatorF::operator()()) to apply the operation to a considered Block
 * 
 * @param qid the reference of the block, see qid_t
 * @param block the Block itself, which cannot be modified
 * @param fid the field on which we operate
 * @param op the ConstOperatorF object containing all the needed data
 */
void ConstCallOpF(const qid_t* qid, GridBlock* block, const Field* fid, ConstOperatorF* op);

//=================================================================================================
/**
 * @brief a Field to Field operator, i.e. which uses the content of one Field to modify another Field
 */
class OperatorF2F {
   public:
    /**
    * @brief Implementation of this virtual function has to be provided by the user as a member function
    * 
    * @warning this function is processed with a multi-thread environment
    * 
    * @param qid the id of the quadrant which corresponds to the current block
    * @param block the current block itself
    * @param fid_src the source field
    * @param fid_trg the target field
    */
    virtual void ApplyOpF2F(const qid_t* qid, GridBlock* block, Field* fid_src, Field* fid_trg) = 0;
    /**
     * @brief call OperatorF2F::ApplyOpF2F() on each block
     */
    virtual void operator()(ForestGrid* grid, Field* field_src, Field* field_trg);
};
/**
 * @brief this function is called by DoOp_() function (through OperatorF2F::operator()()) to apply the operation to a considered Block
 * 
 * @param qid the reference of the block, see qid_t
 * @param block the Block itself, which cannot be modified
 * @param fid_src the source field on which we operate
 * @param fid_trg the traget field on which we operate
 * @param op the OperatorF2F object containing all the needed data
 */
void CallOpF2F(const qid_t* qid, GridBlock* block, const Field* fid_src, Field* fid_trg, OperatorF2F* op);

//=================================================================================================
/**
 * @brief a Field + Field to Field operator, i.e. which uses the content of two Fields to modify another Field
 */
class OperatorFF2F {
   public:
    /**
    * @brief Implementation of this virtual function has to be provided by the user as a member function
    * 
    * @warning this function is processed with a multi-thread environment
    * 
    * @param qid the id of the quadrant which corresponds to the current block
    * @param block the current block itself
    * @param fid_x the source field #1
    * @param fid_y the source field #2
    * @param fid_z the target field
    */
    virtual void ApplyOpFF2F(const qid_t* qid, GridBlock* block, Field* fid_x, Field* fid_y, Field* fid_z) = 0;
    /**
     * @brief call OperatorF2F::ApplyOpF2F() on each block
     */
    virtual void operator()(ForestGrid* grid, Field* field_x, Field* field_y, Field* field_z);
};
/**
 * @brief this function is called by DoOp_() function (through OperatorF2F::operator()()) to apply the operation to a considered Block
 * 
 * @param qid the reference of the block, see qid_t
 * @param block the Block itself, which cannot be modified
 * @param fid_x the source field #1
 * @param fid_y the source field #2
 * @param fid_z the target field
 * @param op the OperatorF2F object containing all the needed data
 */
void CallOpFF2F(const qid_t* qid, GridBlock* block, Field* fid_x, Field* fid_y, Field* fid_z, OperatorF2F* op);

// //=================================================================================================
// /**
//  * @brief defines which mode is triggered in a OperatorDeriv
//  */
// typedef enum mode_deriv_t {
//     M_DERIV_INNER,
//     M_DERIV_OUTER
// } mode_deriv_t;
// /**
//  * @brief a Derivation operator, i.e. which derivates one Field and store the result in a second one
//  *
//  * For this operator, we need to divide the operation into two steps: the inner step, which do not rely on the
//  * ghost points and the outer step, which do rely on ghost points
//  */
// class OperatorDeriv {
//     sid_t src_ida_;  //!< dimension available in the fid_src
//    public:
//     /**
//     * @brief Implementation of this virtual function has to be provided by the user as a member function
//     *
//     * @param qid the id of the quadrant which corresponds to the current block
//     * @param block the current block itself
//     * @param fid_src the source field
//     * @param fid_trg the target field
//     */
//     virtual void ApplyOpDerivInner(const qid_t* qid, GridBlock* block, const Field* fid_src, Field* fid_trg) = 0;
//     virtual void ApplyOpDerivOuter(const qid_t* qid, GridBlock* block, const Field* fid_src, Field* fid_trg) = 0;
//     // defines the operator
//     virtual void operator()(ForestGrid* grid, Field* field_src, Field* field_trg, const sid_t src_ida, const mode_deriv_t mode);
// };
// void CallOpDerivInner(const qid_t* qid, GridBlock* block, const Field* fid_src, Field* fid_trg, OperatorF2F* op);
// void CallOpDerivOuter(const qid_t* qid, GridBlock* block, const Field* fid_src, Field* fid_trg, OperatorF2F* op);

//=================================================================================================
/**
 * @brief a Constant Field and Field operator, i.e. which uses the content of two Fields without modifying it
 * 
 */
class ConstOperatorFF {
   public:
    /**
     * @brief Implementation of this virtual function has to be provided by the user as a member function
     * 
     * @warning this function is processed with a multi-thread environment
     * 
     * @param qid the id of the quadrant which corresponds to the current block
     * @param block the current block itself
     * @param fid_1 a first field
     * @param fid_2 a second field
     */
    virtual void ApplyConstOpFF(const qid_t* qid, GridBlock* block, const Field* fid_1, const Field* fid_2) = 0;
    /**
     * @brief call ConstOperatorFF::ApplyConstOpFF() on each block
     */
    virtual void operator()(ForestGrid* grid, Field* fid_1, Field* fid_2);
};

/**
 * @brief this function is called by DoOp_() function (through ConstOperatorFF::operator()()) to apply the operation to a considered Block
 * 
 * @param qid the reference of the block, see qid_t
 * @param block the Block itself
 * @param fid_1 the first field on which we operate
 * @param fid_2 the second field on which we operate
 * @param op the operator
 */
void ConstCallOpFF(const qid_t* qid, GridBlock* block, const Field* fid_1, const Field* fid_2, ConstOperatorFF* op);

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
    p8est_t*      forest    = grid->forest();
    p8est_mesh_t* mesh      = grid->mesh();
    const iblock_t   nqlocal   = p4est_NumQuadOnLevel(mesh,lvl);

#pragma omp parallel for firstprivate(data)
    for (iblock_t lid = 0; lid < nqlocal; lid++) {
        // get the corresponding id of the quadrant
        const iblock_t bid = p4est_GetQuadIdOnLevel(mesh,lvl,lid);
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

// //=================================================================================================
// /**
//  * @brief wraps the function DoOp_F_() for an operator (`op_t<T, Field*>`, see @ref op_t), which performs the actual implementation
//  *
//  * We first call the operation on the blocks, then set the ghosts as changed
//  *
//  * @tparam T the type of data which is given as context to the operator
//  * @param op the operator itself, see definition of op_t
//  * @param grid the grid on which one want to iterate
//  * @param field the field on which the operation has to be done
//  * @param data the context, i.e. user custom data to forward to the operator
//  */
// template <typename T>
// void DoOp(const op_t<T, Field*> op, ForestGrid* grid, Field* field, T data) {
//     // execute the operator on the blocks
//     DoOp_F_<op_t<T, Field*>, T, Field*>(op, grid, field, data);
//     // set the ghost as changed
//     m_verb("setting the ghosts of %s to false", field->name().c_str());
//     field->ghost_status(false);
// }

// /**
//  * @brief wraps the function DoOp_F_() for a constant operator (`const op_t<T, const Field*>`, see @ref op_t), which performs the actual implementation
//  *
//  *
//  * @tparam T the type of data which is given as context to the operator
//  * @param op the operator itself, see definition of op_t
//  * @param grid the grid on which one want to iterate
//  * @param field the field on which the operation has to be done
//  * @param data the context, i.e. user custom data to forward to the operator
//  */
// template <typename T>
// void DoOp(const op_t<T, const Field*> op, ForestGrid* grid, Field* field, T data) {
//     // execute the operator on the blocks
//     DoOp_F_<op_t<T, const Field*>, T, const Field*>(op, grid, field, data);
//     // no change of the ghost is needed
// }

// /**
//  * @brief wraps the function DoOp_F_() for a field to field operator (`op_t<T, const Field*, Field*>`, see @ref op_t), which performs the actual implementation
//  *
//  * We first call the operation on the blocks, then set the ghosts as changed
//  *
//  * @tparam T the type of data which is given as context to the operator
//  * @param op the operator itself, see definition of op_t
//  * @param grid the grid on which one want to iterate
//  * @param field the field on which the operation has to be done
//  * @param data the context, i.e. user custom data to forward to the operator
//  */
// template <typename T>
// void DoOp(const op_t<T, const Field*, Field*> op, ForestGrid* grid, Field* field_src, Field* field_trg, T data) {
//     // execute the operation on the blocks
//     DoOp_F_<op_t<T, const Field*, Field*>, T, const Field*, Field*>(op, grid, field_src, field_trg, data);
//     // set the field_trg ghosts as changed
//     m_verb("setting the ghosts of %s to false", field_trg->name().c_str());
//     field_trg->ghost_status(false);
// }

#endif  // SRC_OPERATOR_HPP_
