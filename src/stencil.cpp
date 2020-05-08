#include "stencil.hpp"

/**
 * @brief Construct a new Stencil by remembering the grid on which we operate.
 * 
 * @note: this is not needed but helps to defines another operator() function that does not conflict with the OperatorF2F::operator()() function
 * 
 * @param grid 
 */
Stencil::Stencil(Grid* grid) {
    grid_ = grid;
}

/**
 * @brief applies the stencil on the current dimension, ida_, either the inner or the outer given the inner_ value 
 * 
 * @param qid 
 * @param block 
 * @param fid_src 
 * @param fid_trg 
 */
void Stencil::ApplyOpF2F(const qid_t* qid, GridBlock* block, const Field* fid_src, Field* fid_trg) {
    //-------------------------------------------------------------------------
    if (inner_) {
        ApplyOpDerivInner(qid, block, fid_src, fid_trg);
    } else {
        ApplyOpDerivOuter(qid, block, fid_src, fid_trg);
    }
    //-------------------------------------------------------------------------
}

/**
 * @brief apply the stencil on the field_src and store the result in field_trg
 * 
 * This functions implements the dimension ghost computation, inner stencil and outer stencil computation overlapping.
 * The computation is done dimension by dimension on the exchanged dimension of the field_src.
 * At the end, the field_src ghost status is changed to true and the target field contains the result.
 * 
 * @param field_src 
 * @param field_trg 
 */
void Stencil::operator()(Field* field_src, Field* field_trg) {
    m_begin;
    m_assert(field_src != nullptr, "the source field cannot be null");
    //-------------------------------------------------------------------------
    // start the send in the first dimension
    grid_->GhostPullSend(field_src, 0);
    // start the inner operation on the first dimension
    if (field_trg != nullptr) {
        if (grid_->profiler() != nullptr) {
            grid_->profiler()->Start("stencil_inner");
        }
        ida_   = 0;
        inner_ = true;
        OperatorF2F::operator()(grid_, field_src, field_trg);

        if (grid_->profiler() != nullptr) {
            grid_->profiler()->Stop("stencil_inner");
        }
    }
    for (int ida = 1; ida < field_src->lda(); ida++) {
        
        // receive the previous dimension
        grid_->GhostPullRecv(field_src, ida - 1);
        // start the send for the next dimension
        grid_->GhostPullSend(field_src, ida);
        // fill the ghost values of the just-received information
        grid_->GhostPullFill(field_src, ida - 1);
        // do the outer operation if needed with the newly computed ghosts and already do the inner operation for the next dimension
        if (field_trg != nullptr) {
            if (grid_->profiler() != nullptr) {
                grid_->profiler()->Start("stencil_outer");
            }
            // outer operation on the just received dim
            ida_   = ida - 1;
            inner_ = false;
            OperatorF2F::operator()(grid_, field_src, field_trg);
            // new operation on the now received dimension
            if (grid_->profiler() != nullptr) {
                grid_->profiler()->Stop("stencil_outer");
                grid_->profiler()->Start("stencil_inner");
            }
            ida_   = ida;
            inner_ = true;
            OperatorF2F::operator()(grid_, field_src, field_trg);
            if (grid_->profiler() != nullptr) {
                grid_->profiler()->Stop("stencil_inner");
            }
        }
    }
    grid_->GhostPullRecv(field_src, field_src->lda() - 1);
    grid_->GhostPullFill(field_src, field_src->lda() - 1);
    // start the inner operation on the first dimension
    if (field_trg != nullptr) {
        if (grid_->profiler() != nullptr) {
            grid_->profiler()->Start("stencil_outer");
        }
        ida_   = field_src->lda() - 1;
        inner_ = false;
        OperatorF2F::operator()(grid_, field_src, field_trg);
        if (grid_->profiler() != nullptr) {
            grid_->profiler()->Stop("stencil_outer");
        }
    }
    // set that everything is ready for the field
    field_src->ghost_status(true);
    if (field_trg != nullptr) {
        field_trg->ghost_status(false);
    }
    //-------------------------------------------------------------------------
    m_end;
}
