#ifndef SRC_JACOBI_HPP_
#define SRC_JACOBI_HPP_

#include "murphy.hpp"


template <sid_t length>
class Jacobi : public IterativeSolver {
   protected:
    real_t alpha_ = 1.0;
    real_t coef_lapla_[length];  //!< coefficients of the laplacian to apply
   public:
    Jacobi(real_t alpha) {
        // store the relaxation param
        alpha_ = alpha;
        // get the stencil
        if (length == 3) {
            coef_lapla_[0] = +1.0;
            coef_lapla_[1] = -2.0;
            coef_lapla_[2] = +1.0;
        } else if (length == 5) {
            coef_lapla_[0] = -1.0 / 12.0;
            coef_lapla_[1] = +4.0 / 3.0;
            coef_lapla_[2] = -5.0 / 2.0;
            coef_lapla_[3] = +4.0 / 3.0;
            coef_lapla_[4] = -1.0 / 12.0;
        } else {
            m_assert(false, "not coded yet");
        }
    }

    void IterativeSolverPrep(const qid_t* qid, GridBlock* block, Field* fid_sol, Field* fid_rhs, Field* fid_tmp) override;
    void IterativeSolverInner(const qid_t* qid, GridBlock* block, Field* fid_sol, Field* fid_rhs, Field* fid_tmp) override;
    void IterativeSolverOuter(const qid_t* qid, GridBlock* block, Field* fid_sol, Field* fid_rhs, Field* fid_tmp) override;
};

#endif  // SRC_JACOBI_HPP_

#include "jacobi.inc"
