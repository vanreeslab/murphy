#include "stencil.hpp"


Stencil::Stencil(Grid* grid){
    grid_ = grid;
}

void Stencil::ApplyOpF2F(const qid_t* qid, GridBlock* block, const Field* fid_src, Field* fid_trg) {
    //-------------------------------------------------------------------------
    if (inner_) {
        ApplyOpDerivInner(qid, block, fid_src, fid_trg);
    } else {
        ApplyOpDerivOuter(qid, block, fid_src, fid_trg);
    }
    //-------------------------------------------------------------------------
}

void Stencil::operator()(Field* field_src, Field* field_trg) {
    m_begin;
    m_assert(field_src != nullptr, "the source field cannot be null");
    //-------------------------------------------------------------------------
    // start the send in the first dimension
    grid_->GhostPullSend(field_src, 0);
    // start the inner operation on the first dimension
    if (field_trg != nullptr) {
        ida_   = 0;
        inner_ = true;
        OperatorF2F::operator()(grid_, field_src, field_trg);
    }
    for (int ida = 1; ida < field_src->lda(); ida++) {
        // receive the previous dimension
        grid_->GhostPullRecv(field_src, ida - 1);
        // start the send for the next dimension
        grid_->GhostPullSend(field_src, ida);
        // fill the ghost values of the just-received information
        grid_->GhostPullFill(field_src, ida - 1);
        // do the outer operation if needed with the newly computed ghosts
        if (field_trg != nullptr) {
            ida_   = ida - 1;
            inner_ = false;
            OperatorF2F::operator()(grid_, field_src, field_trg);
        }
    }
    grid_->GhostPullRecv(field_src, field_src->lda() - 1);
    // start the inner operation on the first dimension
    if (field_trg != nullptr) {
        ida_   = field_src->lda() - 1;
        inner_ = false;
        OperatorF2F::operator()(grid_, field_src, field_trg);
    }
    // set that everything is ready for the field
    field_src->ghost_status(true);
    if (field_trg != nullptr) {
        field_trg->ghost_status(false);
    }
    //-------------------------------------------------------------------------
    m_end;
}