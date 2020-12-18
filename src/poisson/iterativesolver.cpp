#include "iterativesolver.hpp"

#include "doop.hpp"

void CallIterativeSolverPrep(const qid_t* qid, GridBlock* block, Field* fid_sol, Field* fid_rhs, Field* fid_tmp, IterativeSolver* solver) {
    solver->IterativeSolverPrep(qid, block, fid_sol, fid_rhs, fid_tmp);
}
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
    m_assert(fid_sol->lda() == fid_rhs->lda(),"the solution and the rhs must have the same lda");
#ifndef NDEBUG
    if (fid_tmp != nullptr) {
        m_assert(fid_sol->lda() == fid_tmp->lda(), "the solution and the tmpp must have the same lda");
    }
#endif
    m_assert(false,"need to change this to adapt to the new ghosting");
    //-------------------------------------------------------------------------
    // determine which field will be used to contain the value
    // if a fid_tmp is given, we copy the current state to him and use it for the ghosts (ex Jacobi), if not, we use the sol field (ex: Gauss Seidel)
    // Field* field_to_ghost;
    // bctype_t* tmp_bctype[6];
    // if(fid_tmp != nullptr){
    //     // store the current tmp boundary condition and impose the solution ones
    //     for(iface_t iface=0; iface<6; iface++){
    //         tmp_bctype[iface] = fid_tmp->bctype(iface);
    //         fid_tmp->bctype(fid_sol->bctype(iface),iface);
    //     }
    //     field_to_ghost = fid_tmp;
    //     ida_ = 0;
    //     DoOp_F_<op_t<IterativeSolver*, Field*, Field*, Field*>, IterativeSolver*, Field*, Field*, Field*>(CallIterativeSolverPrep, grid, fid_sol, fid_rhs, fid_tmp, this);
    // }
    // else{
    //     field_to_ghost = fid_sol;
    // }
    // // m_log("doing the ghosts on %s",field_to_ghost->name().c_str());
    
    // // start the send in the first dimension
    // grid->GhostPullSend(field_to_ghost, 0);

    // // start the inner operation on the first dimension
    // m_profStart(grid->profiler(),"stencil_inner");
    // ida_   = 0;
    // DoOp_F_<op_t<IterativeSolver*, Field*, Field*, Field*>, IterativeSolver*, Field*, Field*, Field*>(CallIterativeSolverInner, grid, fid_sol, fid_rhs, fid_tmp, this);
    // m_profStop(grid->profiler(),"stencil_inner");

    // for (int ida = 1; ida < field_to_ghost->lda(); ida++) {
    //     // prepare the work for the next dimension, to get the send ready
    //     if (fid_tmp != nullptr) {
    //         ida_ = ida;
    //         DoOp_F_<op_t<IterativeSolver*, Field*, Field*, Field*>, IterativeSolver*, Field*, Field*, Field*>(CallIterativeSolverPrep, grid, fid_sol, fid_rhs, fid_tmp, this);
    //     }
    //     // receive the previous dimension
    //     grid->GhostPullRecv(field_to_ghost, ida - 1);
    //     // start the send for the next dimension, already prepared
    //     grid->GhostPullSend(field_to_ghost, ida);
    //     // fill the ghost values of the just-received information for the previous dimension
    //     grid->GhostPullFill(field_to_ghost, ida - 1);
    //     // do the outer operation on the previous dimension
    //     m_profStart(grid->profiler(), "stencil_outer");
    //     ida_ = ida - 1;
    //     DoOp_F_<op_t<IterativeSolver*, Field*, Field*, Field*>, IterativeSolver*, Field*, Field*, Field*>(CallIterativeSolverOuter, grid, fid_sol, fid_rhs, fid_tmp, this);
    //     m_profStop(grid->profiler(), "stencil_outer");
    //     // the previous dimension is now over, we can start the new dimension
    //     // do the inner operation on the new dimension
    //     m_profStart(grid->profiler(), "stencil_inner");
    //     ida_ = ida;
    //     DoOp_F_<op_t<IterativeSolver*, Field*, Field*, Field*>, IterativeSolver*, Field*, Field*, Field*>(CallIterativeSolverInner, grid, fid_sol, fid_rhs, fid_tmp, this);
    //     m_profStop(grid->profiler(), "stencil_inner");
    // }
    // grid->GhostPullRecv(field_to_ghost, field_to_ghost->lda() - 1);
    // grid->GhostPullFill(field_to_ghost, field_to_ghost->lda() - 1);
    // // start the inner operation on the first dimension

    // if (grid->profiler() != nullptr) {
    //     grid->profiler->Start("stencil_outer");
    // }
    // ida_ = fid_sol->lda() - 1;
    // DoOp_F_<op_t<IterativeSolver*, Field*, Field*, Field*>, IterativeSolver*, Field*, Field*, Field*>(CallIterativeSolverOuter, grid, fid_sol, fid_rhs, fid_tmp, this);
    // if (grid->profiler() != nullptr) {
    //     grid->profiler->Stop("stencil_outer");
    // }
    // // we changed the solution field
    // fid_sol->ghost_status(false);
    // if (fid_tmp != nullptr) {
    //     // reset the correct boundary conditions
    //     for(iface_t iface=0; iface<6; iface++){
    //         fid_tmp->bctype(tmp_bctype[iface],iface);
    //     }
    //     // set the ghosts to changed
    //     fid_tmp->ghost_status(false);
    // }
    //-------------------------------------------------------------------------
    m_end;
}
