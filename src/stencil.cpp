#include "stencil.hpp"

#include "doop.hpp"

/**
 * @brief Construct a new Stencil with no profiler attached
 * 
 * The Stencil only operates on the inner block and do not fill the ghost points
 * 
 */
Stencil::Stencil() : BlockOperator(nullptr, nullptr) {};

Stencil::Stencil(Prof* profiler) : BlockOperator(nullptr, profiler) {
    m_begin;
    //-------------------------------------------------------------------------
    // m_profCreate(profiler,"stencil_outer");
    // m_profCreate(profiler,"stencil_inner");
    //-------------------------------------------------------------------------
    m_end;
};

/**
 * @brief apply the stencil on the field_src and store the result in field_trg
 * 
 * This functions implements the ghost value update, inner stencil and outer stencil computation overlapping, dimension by dimension.
 * At the end, the field_src ghost status is changed to true and the target field contains the result (with wrong ghost stauts).
 * 
 * @param field_src 
 * @param field_trg 
 * @param grid
 */
void Stencil::operator()(Grid* grid, Field* field_src, Field* field_trg) {
    m_begin;
    m_assert(grid != nullptr, "the grid cannot be null");
    m_assert(field_src != nullptr, "the source field cannot be null");
    m_assert(field_trg != nullptr, "the source field cannot be null");
    m_assert(grid->is_mesh_valid(), "we need the mesh and the ghost to do something here");
    //-------------------------------------------------------------------------
    m_profStart(grid->profiler(),"stencil");
    // init the prof if not already done
    for (lda_t ida = 0; ida < field_src->lda(); ++ida) {
        // start the send of the coarse
        m_profStart(grid->profiler(),"ghost");
        grid->GhostPull_Post(field_src, ida);
        m_profStop(grid->profiler(),"ghost");

        // compute the stencil on the outer side
        m_profStart(grid->profiler(), "inner");
        if (field_trg != nullptr) {
            ida_ = ida;
            DoOpMesh(this,&Stencil::ApplyStencilInner,grid,field_src, field_trg);
        }
        m_profStop(grid->profiler(), "inner");

        // get the coarse representation back
        m_profStart(grid->profiler(),"ghost");
        grid->GhostPull_Wait(field_src, ida);
        m_profStop(grid->profiler(),"ghost");

        // inner operation on the now received dimension
        ida_ = ida;
        m_profStart(grid->profiler(), "outer");
        DoOpMesh(this, &Stencil::ApplyStencilOuter, grid, field_src, field_trg);
        m_profStop(grid->profiler(), "outer");
    }
    m_profStop(grid->profiler(),"stencil");
    // update the ghost status
    field_src->ghost_status(true);
    field_trg->ghost_status(false);
    //-------------------------------------------------------------------------
    m_end;
}

// void Stencil::operator()(Grid* grid, Field* field_src, Field* field_trg) {
//     m_begin;
//     m_assert(field_src != nullptr, "the source field cannot be null");
//     //-------------------------------------------------------------------------
//     // start the send in the first dimension
//     grid->GhostPullSend(field_src, 0);
//     // start the inner operation on the first dimension
//     if (field_trg != nullptr) {
//         ida_ = 0;
//         m_profStart(grid->profiler(), "stencil_inner");
//         DoOp_F_<op_t<Stencil*, Field*, Field*>, Stencil*, Field*, Field*>(CallStencilOpInner, grid, field_src, field_trg, this);
//         m_profStop(grid->profiler(), "stencil_inner");
//     }
//     for (int ida = 1; ida < field_src->lda(); ida++) {
//         // receive the previous dimension
//         grid->GhostPullRecv(field_src, ida - 1);
//         // start the send for the next dimension
//         grid->GhostPullSend(field_src, ida);
//         // fill the ghost values of the just-received information
//         grid->GhostPullFill(field_src, ida - 1);
//         // do the outer operation if needed with the newly computed ghosts and already do the inner operation for the next dimension
//         if (field_trg != nullptr) {
//             // outer operation on the just received dim
//             ida_ = ida - 1;
//             m_profStart(grid->profiler(), "stencil_outer");
//             DoOp_F_<op_t<Stencil*, Field*, Field*>, Stencil*, Field*, Field*>(CallStencilOpOuter, grid, field_src, field_trg, this);
//             m_profStop(grid->profiler(), "stencil_outer");
//             // inner operation on the now received dimension
//             ida_ = ida;
//             m_profStart(grid->profiler(), "stencil_inner");
//             DoOp_F_<op_t<Stencil*, Field*, Field*>, Stencil*, Field*, Field*>(CallStencilOpInner, grid, field_src, field_trg, this);
//             m_profStop(grid->profiler(), "stencil_inner");
//         }
//     }
//     grid->GhostPullRecv(field_src, field_src->lda() - 1);
//     grid->GhostPullFill(field_src, field_src->lda() - 1);
//     // start the inner operation on the first dimension
//     if (field_trg != nullptr) {
//         ida_ = field_src->lda() - 1;
//         m_profStart(grid->profiler(), "stencil_outer");
//         DoOp_F_<op_t<Stencil*, Field*, Field*>, Stencil*, Field*, Field*>(CallStencilOpOuter, grid, field_src, field_trg, this);
//         m_profStop(grid->profiler(), "stencil_outer");
//     }
//     // update the ghost status
//     field_src->ghost_status(true);
//     if (field_trg != nullptr) {
//         field_trg->ghost_status(false);
//     }
//     //-------------------------------------------------------------------------
//     m_end;
// }
