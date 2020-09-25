#ifndef SRC_OPERATOR_HPP_
#define SRC_OPERATOR_HPP_

#include <limits>

#include "doop.hpp"
#include "field.hpp"
#include "forestgrid.hpp"
#include "gridblock.hpp"
#include "murphy.hpp"


using std::nullptr_t;
using std::numeric_limits;

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
// /**
//  * @brief this function is called by DoOp_() function (through OperatorF::operator()()) to apply the operation to a considered Block
//  * 
//  * @param qid the reference of the block, see qid_t
//  * @param block the Block itself
//  * @param fid the field on which we operate
//  * @param op the OperatorF object containing all the needed data
//  */
// static void CallOpF(const qid_t* qid, GridBlock* block, Field* fid, OperatorF* op);

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
    virtual void operator()(ForestGrid* grid, const Field* field);
};
// /**
//  * @brief this function is called by DoOp_() function (through ConstOperatorF::operator()()) to apply the operation to a considered Block
//  * 
//  * @param qid the reference of the block, see qid_t
//  * @param block the Block itself, which cannot be modified
//  * @param fid the field on which we operate
//  * @param op the ConstOperatorF object containing all the needed data
//  */
// static void ConstCallOpF(const qid_t* qid, GridBlock* block, const Field* fid, ConstOperatorF* op);

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
// /**
//  * @brief this function is called by DoOp_() function (through OperatorF2F::operator()()) to apply the operation to a considered Block
//  * 
//  * @param qid the reference of the block, see qid_t
//  * @param block the Block itself, which cannot be modified
//  * @param fid_src the source field on which we operate
//  * @param fid_trg the traget field on which we operate
//  * @param op the OperatorF2F object containing all the needed data
//  */
// static void CallOpF2F(const qid_t* qid, GridBlock* block, const Field* fid_src, Field* fid_trg, OperatorF2F* op);

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
// /**
//  * @brief this function is called by DoOp_() function (through OperatorF2F::operator()()) to apply the operation to a considered Block
//  * 
//  * @param qid the reference of the block, see qid_t
//  * @param block the Block itself, which cannot be modified
//  * @param fid_x the source field #1
//  * @param fid_y the source field #2
//  * @param fid_z the target field
//  * @param op the OperatorF2F object containing all the needed data
//  */
// static void CallOpFF2F(const qid_t* qid, GridBlock* block, Field* fid_x, Field* fid_y, Field* fid_z, OperatorF2F* op);

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
    virtual void operator()(ForestGrid* grid, const Field* fid_1, const Field* fid_2);
};

// /**
//  * @brief this function is called by DoOp_() function (through ConstOperatorFF::operator()()) to apply the operation to a considered Block
//  * 
//  * @param qid the reference of the block, see qid_t
//  * @param block the Block itself
//  * @param fid_1 the first field on which we operate
//  * @param fid_2 the second field on which we operate
//  * @param op the operator
//  */
// static void ConstCallOpFF(const qid_t* qid, GridBlock* block, const Field* fid_1, const Field* fid_2, ConstOperatorFF* op);

#endif  // SRC_OPERATOR_HPP_
