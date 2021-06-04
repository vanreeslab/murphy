#ifndef SRC_OPERATOR_STENCIL_HPP_
#define SRC_OPERATOR_STENCIL_HPP_

#include "core/macros.hpp"
#include "core/types.hpp"
#include "grid/grid.hpp"
#include "operator/blockoperator.hpp"

/**
 * @brief defines a stencil application, that will compute the ghosts on a source field and fill the result in a target field
 * 
 * @warning Because of the overlap with the stencil computation and the ghost exchange for the source field, the application of the stencil
 * is done one dimension of the source field at a time!
 * As an example, if a cross product is needed, we use the current source field dimension, given by @ref ida_, to fill the two other dimensions
 * (given by `(ida_ + 1)%3` and `(ida_ + 2)%3`
 * 
 */
class Stencil : public BlockOperator {
   protected:
    sid_t ida_ = 0;  //!< current source dimension

   public:
    // default void constructor
    explicit Stencil();

    /**
     * @brief returns the number of ghost point needed by the stencil
     * 
     * @return lid_t 
     */
    virtual lid_t NGhost() const = 0;

    /**
     * @brief execute the whole stencil, computation on every block, including the ghost value computation, the inner and outer computation using overlapping
     * between the ghost exchange and the stencil computation.
     */
    void operator()(m_ptr<const Grid> grid, m_ptr<Field> field_src, m_ptr<Field> field_trg);

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
    virtual void DoMagic(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, const bool is_outer, m_ptr<const Field> fid_src, m_ptr<Field> fid_trg) const = 0;
};

#endif  // SRC_OPERATOR_STENCIL_HPP_
