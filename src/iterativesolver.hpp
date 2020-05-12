#ifndef SRC_ITERATIVESOLVER_HPP_
#define SRC_ITERATIVESOLVER_HPP_

#include "grid.hpp"
#include "murphy.hpp"
#include "operator.hpp"

class IterativeSolver : public OperatorF2F {
   protected:
    sid_t ida_ = 0;

   public:
    // //
    // /**
    //  * @brief execute the function ApplyOpDerivInner() or ApplyOpDerivOuter() depending on the inner_ value
    //  */
    virtual void ApplyOpF2F(const qid_t* qid, GridBlock* block, Field* field_sol, Field* field_rhs) = 0;
    /**
     * @brief execute the whole stencil, computation on every block, including the ghost value computation, the inner and outer computation using overlapping
     * between the ghost exchange and the stencil computation.
     * 
     * @note the grid is the last argument to differentiate from the operator() from the @OperatorF2F class
     */
    void operator()(Field* field_src, Field* field_trg, Grid* grid);
};

#endif  // SRC_ITERATIVESOLVER_HPP_