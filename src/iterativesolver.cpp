#include "iterativesolver.hpp"

#include "operator.hpp"

void CallIterativeSolverInner(const qid_t* qid, GridBlock* block, Field* fid_sol, Field* fid_rhs, Field* fid_tmp, IterativeSolver* solver) {
    solver->IterativeSolverInner(qid, block, fid_sol, fid_rhs, fid_tmp);
}
void CallIterativeSolverOuter(const qid_t* qid, GridBlock* block, Field* fid_sol, Field* fid_rhs, Field* fid_tmp, IterativeSolver* solver) {
    solver->IterativeSolverOuter(qid, block, fid_sol, fid_rhs, fid_tmp);
}

void IterativeSolver::operator()(Field* fid_sol, Field* fid_rhs, Field* fid_tmp, Grid* grid) {
    m_begin;
    m_assert(fid_sol != nullptr, "the source field cannot be null");
    m_assert(fid_rhs != nullptr, "the source field cannot be null");
    //-------------------------------------------------------------------------
    // start the send in the first dimension
    grid->GhostPullSend(fid_sol, 0);
    if (grid->profiler() != nullptr) {
        grid->profiler()->Start("stencil_inner");
    }
    ida_   = 0;
    // start the inner operation on the first dimension
    DoOp_F_<op_t<IterativeSolver*, Field*, Field*, Field*>, IterativeSolver*, Field*, Field*, Field*>(CallIterativeSolverInner, grid, fid_sol, fid_rhs, fid_tmp, this);
    if (grid->profiler() != nullptr) {
        grid->profiler()->Stop("stencil_inner");
    }
    for (int ida = 1; ida < fid_sol->lda(); ida++) {
                // receive the previous dimension
        grid->GhostPullRecv(fid_sol, ida - 1);
        // start the send for the next dimension
        grid->GhostPullSend(fid_sol, ida);
        // fill the ghost values of the just-received information
        grid->GhostPullFill(fid_sol, ida - 1);
        // do the outer operation if needed with the newly computed ghosts and already do the inner operation for the next dimension
        if (grid->profiler() != nullptr) {
            grid->profiler()->Start("stencil_outer");
        }
        // outer operation on the just received dim
        ida_   = ida - 1;
        DoOp_F_<op_t<IterativeSolver*, Field*, Field*, Field*>, IterativeSolver*, Field*, Field*, Field*>(CallIterativeSolverOuter ,grid, fid_sol, fid_rhs, fid_tmp, this);
        // new operation on the now received dimension
        if (grid->profiler() != nullptr) {
            grid->profiler()->Stop("stencil_outer");
            grid->profiler()->Start("stencil_inner");
        }
        ida_   = ida;
        DoOp_F_<op_t<IterativeSolver*, Field*, Field*, Field*>, IterativeSolver*, Field*, Field*, Field*>(CallIterativeSolverInner, grid, fid_sol, fid_rhs, fid_tmp, this);
        if (grid->profiler() != nullptr) {
            grid->profiler()->Stop("stencil_inner");
        }
    }
    grid->GhostPullRecv(fid_sol, fid_sol->lda() - 1);
    grid->GhostPullFill(fid_sol, fid_sol->lda() - 1);
    // start the inner operation on the first dimension

    if (grid->profiler() != nullptr) {
        grid->profiler()->Start("stencil_outer");
    }
    ida_ = fid_sol->lda() - 1;
    DoOp_F_<op_t<IterativeSolver*, Field*, Field*, Field*>, IterativeSolver*, Field*, Field*, Field*>(CallIterativeSolverOuter, grid, fid_sol, fid_rhs, fid_tmp, this);
    if (grid->profiler() != nullptr) {
        grid->profiler()->Stop("stencil_outer");
    }
    // set that everything is ready for the field
    fid_sol->ghost_status(false);
    //-------------------------------------------------------------------------
    m_end;
}
