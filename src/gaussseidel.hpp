#ifndef SRC_GAUSSSEIDEL_HPP_
#define SRC_GAUSSSEIDEL_HPP_


#include "murphy.hpp"
#include "laplacian.hpp"

template <sid_t length>
class GaussSeidel : public LaplacianCross<length> {
   protected:
    real_t alpha_ = 1.0;  //!< over-relaxation parameter (typically for Jacobi or GS)
   public:
    GaussSeidel(real_t alpha) : LaplacianCross<length>() {
        // store the relaxation param
        alpha_ = alpha;
    }

   protected:
    virtual void ApplyOpDerivInner(const qid_t* qid, GridBlock* block, const Field* fid_rhs, Field* fid_sol) override;
    virtual void ApplyOpDerivOuter(const qid_t* qid, GridBlock* block, const Field* fid_rhs, Field* fid_sol) override;
};

#endif  // SRC_GAUSSSEIDEL_HPP_

#include "gaussseidel.ipp"