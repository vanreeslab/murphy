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
void Stencil::operator()(m_ptr<const Grid> grid, m_ptr<Field> field_src, m_ptr<Field> field_trg){
    m_begin;
    m_assert(!grid.IsEmpty(), "the grid cannot be null");
    m_assert(!field_src.IsEmpty(), "the source field cannot be null");
    m_assert(!field_trg.IsEmpty(), "the source field cannot be null");
    m_assert(grid->is_mesh_valid(), "we need the mesh and the ghost to do something here");
    m_assert(grid->NGhostFront() >= this->NGhost(), "the wavelet do not provied enough ghost points for the stencil: %d vs %d", grid->NGhostFront(), this->NGhost());
    m_assert(grid->NGhostBack() >= this->NGhost(), "the wavelet do not provied enough ghost points for the stencil: %d vs %d", grid->NGhostBack(), this->NGhost());
    //-------------------------------------------------------------------------
    m_log("ghost check: field <%s> is %s", field_src->name().c_str(), field_src->ghost_status() ? "OK" : "to be computed");
    m_profStart(prof_(), "stencil");
    // init the prof if not already done
    for (lda_t ida = 0; ida < field_src->lda(); ++ida) {
        // start the send of the coarse
        m_profStart(prof_(), "ghost");
        grid->GhostPull_Post(field_src, ida);
        m_profStop(prof_(), "ghost");

        // compute the stencil on the outer side
        m_profStart(prof_(), "inner");
        if (!field_trg.IsEmpty()) {
            ida_ = ida;
            // DoOpMesh(this, &Stencil::ApplyStencilInner, grid, field_src, field_trg);
            DoOpMesh(this, &Stencil::DoMagic, grid, false, field_src, field_trg);
        }
        m_profStop(prof_(), "inner");

        // get the coarse representation back
        m_profStart(prof_(), "ghost");
        grid->GhostPull_Wait(field_src, ida);
        m_profStop(prof_(), "ghost");

        // outer operation on the now received dimension
        ida_ = ida;
        m_profStart(prof_(), "outer");
        // DoOpMesh(this, &Stencil::ApplyStencilOuter, grid, field_src, field_trg);
        DoOpMesh(this, &Stencil::DoMagic, grid, true, field_src, field_trg);
        m_profStop(prof_(), "outer");
    }
    m_profStop(prof_(), "stencil");
    // update the ghost status
    field_src->ghost_status(true);
    field_trg->ghost_status(false);
    //-------------------------------------------------------------------------
    m_end;
}