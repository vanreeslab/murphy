#include "operator.hpp"

//=================================================================================================
void CallOpF(const qid_t* qid, Block* block, Field* fid, OperatorF* op) {
    op->apply(qid, block, fid);
}
void OperatorF::DoOp(Grid* grid, Field* field) {
    // execute the operator on the blocks
    DoOp_<op_t<OperatorF*, Field*>, OperatorF*, Field*>(CallOpF, grid, field, this);
    // set the ghost as changed
    m_verb("setting the ghosts of %s to false", field->name().c_str());
    field->SetGhostStatus(false);
}

//=================================================================================================
void ConstCallOpF(const qid_t* qid, Block* block, const Field* fid, ConstOperatorF* op) {
    op->apply(qid, block, fid);
}
void ConstOperatorF::DoOp(Grid* grid, Field* field) {
    // execute the operation on the blocks
    DoOp_<op_t<ConstOperatorF*, const Field*>, ConstOperatorF*, const Field*>(ConstCallOpF, grid, field, this);
}
//=================================================================================================
void CallOpF2F(const qid_t* qid, Block* block, const Field* fid_src, Field* fid_trg, OperatorF2F* op) {
    op->apply(qid, block, fid_src, fid_trg);
}
void OperatorF2F::DoOp(Grid* grid, const Field* field_src, Field* field_trg) {
    // execute the operation on the blocks
    DoOp_<op_t<OperatorF2F*, const Field*, Field*>, OperatorF2F*, const Field*, Field*>(CallOpF2F, grid, field_src, field_trg, this);
    // set the field_trg ghosts as changed
    m_verb("setting the ghosts of %s to false", field_trg->name().c_str());
    field_trg->SetGhostStatus(false);
}