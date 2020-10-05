// #include "operator.hpp"

// //=================================================================================================
// static void CallOpF(const qid_t* qid, GridBlock* block, Field* fid, OperatorF* op) {
//     op->ApplyOpF(qid, block, fid);
// }
// void OperatorF::operator()(ForestGrid* grid, Field* field) {
//     // execute the operator on the blocks
//     DoOpMesh(&CallOpF, grid, field, this);
//     // set the ghost as changed
//     m_verb("setting the ghosts of %s to false", field->name().c_str());
//     field->ghost_status(false);
// }

// //=================================================================================================
// static void ConstCallOpF(const qid_t* qid, GridBlock* block, const Field* fid, ConstOperatorF* op) {
//     op->ApplyConstOpF(qid, block, fid);
// }
// void ConstOperatorF::operator()(ForestGrid* grid, const Field* field) {
//     // execute the operation on the blocks
//     // DoOp_F_<op_t<ConstOperatorF*, const Field*>, ConstOperatorF*, const Field*>(ConstCallOpF, grid, field, this);
//     DoOpMesh(&ConstCallOpF, grid, field, this);
// }

// //=================================================================================================
// void CallOpF2F(const qid_t* qid, GridBlock* block, Field* fid_src, Field* fid_trg, OperatorF2F* op) {
//     op->ApplyOpF2F(qid, block, fid_src, fid_trg);
// }
// /**
//  * @brief starts the call on every block, i.e. call DoOp_F_ that in turns will call ApplyOpF2F on each block
//  * 
//  * @param grid the grid on which we iterate
//  * @param field_src the source field for the computation, its ghost status is NOT changed
//  * @param field_trg the target field for the computation, its ghost status IS changed
//  */
// void OperatorF2F::operator()(ForestGrid* grid, Field* field_src, Field* field_trg) {
//     // execute the operation on the blocks
//     // DoOp_F_<op_t<OperatorF2F*, Field*, Field*>, OperatorF2F*, Field*, Field*>(CallOpF2F, grid, field_src, field_trg, this);
//     DoOpMesh(&CallOpF2F, grid, field_src, field_trg, this);
//     // set the field_trg ghosts as changed
//     m_verb("setting the ghosts of %s to false", field_trg->name().c_str());
//     field_trg->ghost_status(false);
// }

// //=================================================================================================
// static void CallOpFF2F(const qid_t* qid, GridBlock* block, Field* fid_x, Field* fid_y, Field* fid_z, OperatorFF2F* op) {
//     op->ApplyOpFF2F(qid, block, fid_x, fid_y, fid_z);
// }
// /**
//  * @brief starts the call on every block, i.e. call DoOp_F_ that in turns will call ApplyOpF2F on each block
//  * 
//  * @param grid the grid on which we iterate
//  * @param field_src the source field for the computation, its ghost status is NOT changed
//  * @param field_trg the target field for the computation, its ghost status IS changed
//  */
// void OperatorFF2F::operator()(ForestGrid* grid, Field* field_x, Field* field_y, Field* field_z) {
//     // execute the operation on the blocks
//     // DoOp_F_<op_t<OperatorFF2F*, Field*, Field*, Field*>, OperatorFF2F*, Field*, Field*, Field*>(CallOpFF2F, grid, field_x, field_y, field_z, this);
//     DoOpMesh(&CallOpFF2F, grid, field_x, field_y, field_z, this);
//     // set the field_trg ghosts as changed
//     m_verb("setting the ghosts of %s to false", field_z->name().c_str());
//     field_z->ghost_status(false);
// }

// //=================================================================================================
// static void ConstCallOpFF(const qid_t* qid, GridBlock* block, const Field* fid_1, const Field* fid_2, ConstOperatorFF* op) {
//     op->ApplyConstOpFF(qid, block, fid_1, fid_2);
// }
// void ConstOperatorFF::operator()(ForestGrid* grid, const Field* fid_1, const Field* fid_2) {
//     // execute the operation on the blocks
//     // DoOp_F_<op_t<ConstOperatorFF*, const Field*, const Field*>, ConstOperatorFF*, const Field*, const Field*>(ConstCallOpFF, grid, fid_1, fid_2, this);
//     DoOpMesh(&ConstCallOpFF, grid, fid_1, fid_2, this);
// }