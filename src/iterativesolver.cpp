#include "iterativesolver.hpp"


void IterativeSolver::operator()(Field* field_sol,Field*field_rhs,Grid* grid){
    m_begin;
    m_assert(field_sol != nullptr, "the source field cannot be null");
    m_assert(field_rhs != nullptr, "the source field cannot be null");
    //-------------------------------------------------------------------------
    // start the send in the first dimension
    grid->GhostPullSend(field_sol, 0);
    // start the inner operation on the first dimension
    
        // if (grid->profiler() != nullptr) {
        //     grid->profiler()->Start("stencil_inner");
        // }
        // ida_   = 0;
        // OperatorF2F::operator()(grid, field_sol, field_rhs);
        // if (grid->profiler() != nullptr) {
        //     grid->profiler()->Stop("stencil_inner");
        // }
    // for (int ida = 1; ida < field_sol->lda(); ida++) {
    //             // receive the previous dimension
    //     grid->GhostPullRecv(field_src, ida - 1);
    //     // start the send for the next dimension
    //     grid->GhostPullSend(field_src, ida);
    //     // fill the ghost values of the just-received information
    //     grid->GhostPullFill(field_src, ida - 1);
    //     // do the outer operation if needed with the newly computed ghosts and already do the inner operation for the next dimension
    //     if (field_trg != nullptr) {
    //         if (grid->profiler() != nullptr) {
    //             grid->profiler()->Start("stencil_outer");
    //         }
    //         // outer operation on the just received dim
    //         ida_   = ida - 1;
    //         inner_ = false;
    //         OperatorF2F::operator()(grid, field_src, field_trg);
    //         // new operation on the now received dimension
    //         if (grid->profiler() != nullptr) {
    //             grid->profiler()->Stop("stencil_outer");
    //             grid->profiler()->Start("stencil_inner");
    //         }
    //         ida_   = ida;
    //         inner_ = true;
    //         OperatorF2F::operator()(grid, field_src, field_trg);
    //         if (grid->profiler() != nullptr) {
    //             grid->profiler()->Stop("stencil_inner");
    //         }
    //     }
    // }
    grid->GhostPullRecv(field_sol, field_sol->lda() - 1);
    grid->GhostPullFill(field_sol, field_sol->lda() - 1);
    // start the inner operation on the first dimension
    
        // if (grid->profiler() != nullptr) {
        //     grid->profiler()->Start("stencil_outer");
        // }
        ida_   = field_sol->lda() - 1;
        // inner_ = false;
        OperatorF2F::operator()(grid, field_sol, field_rhs);
    //     if (grid->profiler() != nullptr) {
    //         grid->profiler()->Stop("stencil_outer");
    //     }
    // }
    // set that everything is ready for the field
    // note: the order is REALLY important, especially if this is a gauss-seidel call...
    field_sol->ghost_status(false);
    //-------------------------------------------------------------------------
    m_end;
}