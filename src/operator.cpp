#include "operator.hpp"



//=================================================================================================
void CallOpS(const qid_t* qid, GridBlock* block, Field* fid, OperatorS* op) {
    m_assert(fid == nullptr, "the field has to be NULL for an OperatorS object");
    op->ApplyOperatorS(qid, block);
}
void OperatorS::operator()(ForestGrid* grid) {
    // execute the operator on the blocks, no 
    DoOp_F_<op_t<OperatorS*, Field*>, OperatorS*, Field*>(CallOpS, grid, nullptr, this);
}


//=================================================================================================
void CallOpF(const qid_t* qid, GridBlock* block, Field* fid, OperatorF* op) {
    op->ApplyOperatorF(qid, block, fid);
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
    op->ApplyConstOperatorF(qid, block, fid);
}
void ConstOperatorF::operator()(ForestGrid* grid, Field* field) {
    // execute the operation on the blocks
    DoOp_F_<op_t<ConstOperatorF*, const Field*>, ConstOperatorF*, const Field*>(ConstCallOpF, grid, field, this);
}
//=================================================================================================
void CallOpF2F(const qid_t* qid, GridBlock* block, const Field* fid_src, Field* fid_trg, OperatorF2F* op) {
    op->ApplyOperatorF2F(qid, block, fid_src, fid_trg);
}
void OperatorF2F::operator()(ForestGrid* grid, const Field* field_src, Field* field_trg) {
    // execute the operation on the blocks
    DoOp_F_<op_t<OperatorF2F*, const Field*, Field*>, OperatorF2F*, const Field*, Field*>(CallOpF2F, grid, field_src, field_trg, this);
    // set the field_trg ghosts as changed
    m_verb("setting the ghosts of %s to false", field_trg->name().c_str());
    field_trg->ghost_status(false);
}