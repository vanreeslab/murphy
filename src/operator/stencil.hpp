#ifndef SRC_OPERATOR_STENCIL_HPP_
#define SRC_OPERATOR_STENCIL_HPP_

#include "core/macros.hpp"
#include "core/types.hpp"
#include "grid/grid.hpp"
#include "grid/gridblock.hpp"
#include "operator/blockoperator.hpp"
#include "doop.hpp"

/**
 * @brief defines a stencil application, that will compute the ghosts on a source field and fill the result in a target field
 * 
 * @warning Because of the overlap with the stencil computation and the ghost exchange for the source field, the application of the stencil
 * is done one dimension of the source field at a time!
 * As an example, if a cross product is needed, we use the current source field dimension, given by @ref ida_, to fill the two other dimensions
 * (given by `(ida_ + 1)%3` and `(ida_ + 2)%3`
 * 
 */
template<typename BlockType>
class Stencil : public BlockOperator {
   protected:
    lda_t ida_ = 0;  //!< current source dimension

   public:
    /**
    * @brief Construct a new Stencil with no profiler attached
    * 
    * The Stencil only operates on the inner block and do not fill the ghost points
    * 
    */
    explicit Stencil(): BlockOperator(nullptr){};

    /**
     * @brief apply the stencil on the field_src and store the result in the field_trg
     *
     * This functions implements the ghost value update, inner stencil and outer stencil computation overlapping, dimension by dimension.
     * At the end, the field_src ghost status is changed to true and the target field contains the result (with wrong ghost status).
     *
     */
    void operator()(const Grid* grid, Field* field_src, Field* field_trg){
        m_begin;
        m_assert(!(grid == nullptr), "the grid cannot be null");
        // m_assert(!(field_src == nullptr), "the source field cannot be null");
        // m_assert(!(field_trg == nullptr), "the source field cannot be null"); 
        m_assert(grid->is_mesh_valid(), "we need the mesh and the ghost to do something here");
        m_assert((ghost_len_res_[0] == 0) && (ghost_len_res_[1] == 0), "the ghost_len_res must be 0 0 instead of %d %d", ghost_len_res_[0], ghost_len_res_[1]);
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
            ida_ = ida;
            DoOpMesh(this, &Stencil::DoMagic, grid, false, field_src, field_trg);
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
        field_src->ghost_len(ghost_len_actual);

        // the trg has been overwritten anyway
        if (!(field_trg == nullptr)) {
            field_trg->ghost_len(ghost_len_res_);
        }
        //-------------------------------------------------------------------------
        m_end;
    };

    /**
     * @brief applies the magic of the stencil on the field_src, in the dimension ida_ only! (inner computation or outer computation depending on is_outer)
     * 
     * @warning for the inner computation, no ghost point can be used.
     * @warning only the dimension ida_ of the source field can be used for the outer computations
     * 
     * @param qid the id of the block
     * @param block the id of the block
     * @param is_outer indicate if we want to apply the inner or the outer part of the stencil
     * @param fid_src the source field, only its dimension ida_ should be used
     * @param fid_trg the target field where all dimensions can be filled
     */
    virtual void DoMagic(const qid_t*  qid, BlockType*  block, const bool is_outer, const Field*  fid_src, Field*  fid_trg) const = 0;
};

// Explicit instantiation of the class we will use 
template class Stencil<GridBlock>;

#endif  // SRC_OPERATOR_STENCIL_HPP_
