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
    // default void constructor
    explicit Stencil();

    /**
     * @brief execute the whole stencil, computation on every block, including the ghost value computation, the inner and outer computation using overlapping
     * between the ghost exchange and the stencil computation.
     */
    void operator()(const Grid* grid, Field* field_src, Field* field_trg);

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

#endif  // SRC_OPERATOR_STENCIL_HPP_
