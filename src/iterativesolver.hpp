#ifndef SRC_ITERATIVESOLVER_HPP_
#define SRC_ITERATIVESOLVER_HPP_

#include "grid.hpp"
#include "defs.hpp"

class IterativeSolver{
   protected:
    sid_t ida_   = 0;

   public:
    /**
     * @brief execute the whole stencil, computation on every block, including the ghost value computation, the inner and outer computation using overlapping
     * between the ghost exchange and the stencil computation.
     * 
     * @note the grid is the last argument to differentiate from the operator() from the @OperatorF2F class
     */
    void operator()(Field* field_sol, Field* field_rhs, Field* field_tmp, Grid* grid);

    virtual void IterativeSolverPrep(const qid_t* qid, GridBlock* block, Field* fid_sol, Field* fid_rhs, Field* field_tmp) = 0;
    virtual void IterativeSolverInner(const qid_t* qid, GridBlock* block, Field* fid_sol, Field* fid_rhs, Field* field_tmp)  = 0;
    virtual void IterativeSolverOuter(const qid_t* qid, GridBlock* block, Field* fid_sol, Field* fid_rhs, Field* field_tmp)  = 0;
};

#endif  // SRC_ITERATIVESOLVER_HPP_