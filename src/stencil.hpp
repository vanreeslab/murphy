#ifndef SRC_STENCIL_HPP_
#define SRC_STENCIL_HPP_

#include "blockoperator.hpp"
#include "grid.hpp"
#include "murphy.hpp"

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
    sid_t ida_  = 0;  //!< current source dimension
    
   public:
    // default void constructor
    explicit Stencil();

    /**
     * @brief execute the whole stencil, computation on every block, including the ghost value computation, the inner and outer computation using overlapping
     * between the ghost exchange and the stencil computation.
     */
    void operator()(Grid* grid, Field* field_src, Field* field_trg);

    /**
    * @brief applies the inner part of the stencil on the field_src, in the dimension ida_ only, i.e. the part that doesn't require ghost values
    * 
    * @param qid the id of the block
    * @param block the GridBlock on which we execute
    * @param fid_src the source field, only its dimension ida_ should be used
    * @param fid_trg the target field where all dimensions can be filled
    */
    virtual void ApplyStencilInner(const qid_t* qid, GridBlock* block, Field* fid_src, Field* fid_trg) = 0;

    /**
     * @brief applies the outer part of the stencil on the field_src, in the dimension ida_ only, i.e. the part that does require ghost values
     * 
     * @param qid the ID of the block to fill
     * @param block the considered block
     * @param fid_src the source field, only its dimension ida_ MUST BE used (the others don't have valid ghost point values!!)
     * @param fid_trg the target field, every dimension can be filled
     */
    virtual void ApplyStencilOuter(const qid_t* qid, GridBlock* block, Field* fid_src, Field* fid_trg) = 0;
};

#endif  // SRC_STENCIL_HPP_
