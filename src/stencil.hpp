#ifndef SRC_STENCIL_HPP_
#define SRC_STENCIL_HPP_

#include "grid.hpp"
#include "murphy.hpp"
#include "operator.hpp"

/**
 * @brief defines a stencil application, a particular type of OperatorF2F that needs the ghost values to execute
 * 
 * @warning Because of the overlap with the stencil computation and the ghost exchange for the source field, the application of the stencil
 * is done one dimension of the source field at a time!
 * As an example, if a cross product is needed, we use the current source field dimension, given by @ref ida_, to fill the two other dimensions
 * (given by `(ida_ + 1)%3` and `(ida_ + 2)%3`
 * 
 */
class Stencil : public OperatorF2F {
   protected:
    sid_t ida_   = 0;     //!< current source dimension
    bool  inner_ = true;  //!< application mode: true = inner, false = outer

    Grid* grid_;

   public:
    explicit Stencil(Grid* grid);

    //
    /**
     * @brief execute the function ApplyOpDerivInner() or ApplyOpDerivOuter() depending on the inner_ value
     */
    void ApplyOpF2F(const qid_t* qid, GridBlock* block, const Field* fid_src, Field* fid_trg) override;
    /**
     * @brief execute the whole stencil, computation on every block, including the ghost value computation, the inner and outer computation using overlapping
     * between the ghost exchange and the stencil computation.
     */
    void operator()(Field* field_src, Field* field_trg);

   protected:
    /**
    * @brief applies the inner part of the stencil on the field_src, in the dimension ida_ only, i.e. the part that doesn't require ghost values
    * 
    * @param qid the id of the block
    * @param block the GridBlock on which we execute
    * @param fid_src the source field, only its dimension ida_ should be used
    * @param fid_trg the target field where all dimensions can be filled
    */
    virtual void ApplyOpDerivInner(const qid_t* qid, GridBlock* block, const Field* fid_src, Field* fid_trg) = 0;

    /**
     * @brief applies the outer part of the stencil on the field_src, in the dimension ida_ only, i.e. the part that does require ghost values
     * 
     * @param qid the ID of the block to fill
     * @param block the considered block
     * @param fid_src the source field, only its dimension ida_ MUST BE used (the others don't have valid ghost point values!!)
     * @param fid_trg the target field, every dimension can be filled
     */
    virtual void ApplyOpDerivOuter(const qid_t* qid, GridBlock* block, const Field* fid_src, Field* fid_trg) = 0;
};

#endif  // SRC_STENCIL_HPP_
