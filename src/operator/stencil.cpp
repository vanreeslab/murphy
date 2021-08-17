#include "stencil.hpp"

#include "doop.hpp"

/**
 * @brief Construct a new Stencil with no profiler attached
 * 
 * The Stencil only operates on the inner block and do not fill the ghost points
 * 
 */
Stencil::Stencil() : BlockOperator(nullptr){};

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
void Stencil::operator()(const Grid* grid, Field* field_src, Field* field_trg) {
    m_begin;
    m_assert(!(grid == nullptr), "the grid cannot be null");
    // m_assert(!(field_src == nullptr), "the source field cannot be null");
    // m_assert(!(field_trg == nullptr), "the source field cannot be null");
    m_assert(grid->is_mesh_valid(), "we need the mesh and the ghost to do something here");
    m_assert((ghost_len_res_[0] == 0) && (ghost_len_res_[1] == 0), "the ghost_len_res must be 0 0 instead of %d %d", ghost_len_res_[0], ghost_len_res_[1]);
    // m_assert(grid->NGhostFront() >= this->NGhost(), "the wavelet do not provied enough ghost points for the stencil: %d vs %d", grid->NGhostFront(), this->NGhost());
    // m_assert(grid->NGhostBack() >= this->NGhost(), "the wavelet do not provied enough ghost points for the stencil: %d vs %d", grid->NGhostBack(), this->NGhost());
    //-------------------------------------------------------------------------
    bidx_t ghost_len[2] = {ghost_len_need_[0], ghost_len_need_[1]};
    grid->GhostPull_SetLength(field_src, ghost_len);
    // m_log("ghost check: field <%s> is %s", field_src->name().c_str(), field_src->ghost_status(ghost_len) ? "OK" : "to be computed");

    m_profStart(prof_, "stencil");
    // init the prof if not already done
    for (lda_t ida = 0; ida < field_src->lda(); ++ida) {
        // start the send of the coarse
        m_profStart(prof_, "ghost");
        grid->GhostPull_Post(field_src, ida, ghost_len);
        m_profStop(prof_, "ghost");

        // compute the stencil on the outer side
        m_profStart(prof_, "inner");
        // if (!(field_trg == nullptr)) {
            ida_ = ida;
            DoOpMesh(this, &Stencil::DoMagic, grid, false, field_src, field_trg);
        // }
        m_profStop(prof_, "inner");

        // get the coarse representation back
        m_profStart(prof_, "ghost");
        grid->GhostPull_Wait(field_src, ida, ghost_len);
        m_profStop(prof_, "ghost");

        // outer operation on the now received dimension
        ida_ = ida;
        m_profStart(prof_, "outer");
        DoOpMesh(this, &Stencil::DoMagic, grid, true, field_src, field_trg);
        m_profStop(prof_, "outer");
    }
    m_profStop(prof_, "stencil");

    // update the ghost status
    // we might have not ghosted the field as we already has some up-to-date information
    const bidx_t ghost_len_actual[2] = {m_max(ghost_len[0], field_src->get_ghost_len(0)),
                                        m_max(ghost_len[1], field_src->get_ghost_len(1))};
    field_src.ghost_len(ghost_len_actual);

    // the trg has been overwritten anyway
    field_trg.ghost_len(ghost_len_res_);
    //-------------------------------------------------------------------------
    m_end;
}