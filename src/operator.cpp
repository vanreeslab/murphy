#include "operator.hpp"

//=================================================================================================
void CallOpS(const qid_t* qid, GridBlock* block, Field* fid, OperatorS* op) {
    m_assert(fid == nullptr, "the field has to be NULL for an OperatorS object");
    op->ApplyOpS(qid, block);
}
void OperatorS::operator()(ForestGrid* grid) {
    // execute the operator on the blocks, no 
    DoOp_F_<op_t<OperatorS*, Field*>, OperatorS*, Field*>(CallOpS, grid, nullptr, this);
}

//=================================================================================================
void CallOpF(const qid_t* qid, GridBlock* block, Field* fid, OperatorF* op) {
    op->ApplyOpF(qid, block, fid);
}
void OperatorF::operator()(ForestGrid* grid, Field* field) {
    // execute the operator on the blocks
    DoOp_F_<op_t<OperatorF*, Field*>, OperatorF*, Field*>(CallOpF, grid, field, this);
    // set the ghost as changed
    m_verb("setting the ghosts of %s to false", field->name().c_str());
    field->ghost_status(false);
}

//=================================================================================================
void ConstCallOpF(const qid_t* qid, GridBlock* block, const Field* fid, ConstOperatorF* op) {
    op->ApplyConstOpF(qid, block, fid);
}
void ConstOperatorF::operator()(ForestGrid* grid, Field* field) {
    // execute the operation on the blocks
    DoOp_F_<op_t<ConstOperatorF*, const Field*>, ConstOperatorF*, const Field*>(ConstCallOpF, grid, field, this);
}

//=================================================================================================
void CallOpF2F(const qid_t* qid, GridBlock* block, const Field* fid_src, Field* fid_trg, OperatorF2F* op) {
    op->ApplyOpF2F(qid, block, fid_src, fid_trg);
}
/**
 * @brief starts the call on every block, i.e. call DoOp_F_ that in turns will call ApplyOpF2F on each block
 * 
 * @param grid the grid on which we iterate
 * @param field_src the source field for the computation, its ghost status is NOT changed
 * @param field_trg the target field for the computation, its ghost status IS changed
 */
void OperatorF2F::operator()(ForestGrid* grid, Field* field_src, Field* field_trg) {
    // execute the operation on the blocks
    DoOp_F_<op_t<OperatorF2F*, const Field*, Field*>, OperatorF2F*, const Field*, Field*>(CallOpF2F, grid, field_src, field_trg, this);
    // set the field_trg ghosts as changed
    m_verb("setting the ghosts of %s to false", field_trg->name().c_str());
    field_trg->ghost_status(false);
}

//=================================================================================================
void CallOpFF2F(const qid_t* qid, GridBlock* block, const Field* fid_x, const Field* fid_y, Field* fid_z, OperatorFF2F* op) {
    op->ApplyOpFF2F(qid, block, fid_x, fid_y, fid_z);
}
/**
 * @brief starts the call on every block, i.e. call DoOp_F_ that in turns will call ApplyOpF2F on each block
 * 
 * @param grid the grid on which we iterate
 * @param field_src the source field for the computation, its ghost status is NOT changed
 * @param field_trg the target field for the computation, its ghost status IS changed
 */
void OperatorFF2F::operator()(ForestGrid* grid, Field* field_x, Field* field_y, Field* field_z) {
    // execute the operation on the blocks
    DoOp_F_<op_t<OperatorFF2F*, const Field*,const Field*, Field*>, OperatorFF2F*, const Field*,const Field*, Field*>(CallOpFF2F, grid, field_x, field_y, field_z, this);
    // set the field_trg ghosts as changed
    m_verb("setting the ghosts of %s to false", field_z->name().c_str());
    field_z->ghost_status(false);
}

// //=================================================================================================
// void CallOpDerivInner(const qid_t* qid, GridBlock* block, const Field* fid_src, Field* fid_trg, OperatorDeriv* op) {
//     op->ApplyOpDerivInner(qid, block, fid_src, fid_trg);
// }
// void CallOpDerivOuter(const qid_t* qid, GridBlock* block, const Field* fid_src, Field* fid_trg, OperatorDeriv* op) {
//     op->ApplyOpDerivOuter(qid, block, fid_src, fid_trg);
// }
// /**
//  * @brief starts the call on every block, i.e. call DoOp_F_ that in turns will call ApplyOpDerivInner or ApplyOpDerivOuter on each block
//  * 
//  * @param grid the grid on which we iterate
//  * @param field_src the source field for the computation, its ghost status is NOT changed
//  * @param field_trg the target field for the computation, its ghost status IS changed
//  * @param src_ida indicate the available dimension in the source field
//  * @param mode the mode to call on the deriv operator
//  */
// void OperatorDeriv::operator()(ForestGrid* grid, Field* field_src, Field* field_trg,const sid_t src_ida, const mode_deriv_t mode) {
//     // register the dimension available in the source field
//     src_ida_ = src_ida;
//     // execute the operation on the blocks
//     if (mode == M_DERIV_INNER) {
//         DoOp_F_<op_t<OperatorDeriv*, const Field*, Field*>, OperatorDeriv*, const Field*, Field*>(CallOpDerivInner, grid, field_src, field_trg, this);
//     } else if (mode == M_DERIV_OUTER) {
//         DoOp_F_<op_t<OperatorDeriv*, const Field*, Field*>, OperatorDeriv*, const Field*, Field*>(CallOpDerivOuter, grid, field_src, field_trg, this);
//     }
//     // set the field_trg ghosts as changed
//     m_verb("setting the ghosts of %s to false", field_trg->name().c_str());
//     field_trg->ghost_status(false);
// }

//=================================================================================================
void ConstCallOpFF(const qid_t* qid, GridBlock* block, const Field* fid_1, const Field* fid_2, ConstOperatorFF* op) {
    op->ApplyConstOpFF(qid, block, fid_1, fid_2);
}
void ConstOperatorFF::operator()(ForestGrid* grid, Field* fid_1, Field* fid_2) {
    // execute the operation on the blocks
    DoOp_F_<op_t<ConstOperatorFF*, const Field*, const Field*>, ConstOperatorFF*, const Field*, const Field*>(ConstCallOpFF, grid, fid_1, fid_2, this);
}