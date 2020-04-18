#include "setvalues.hpp"

SetGaussian::SetGaussian(real_t sigma, real_t center[3]){
    m_begin;
    //-------------------------------------------------------------------------
    sigma_ = sigma;
    for(int id=0; id<3; id++){
        center_[id] = center[id];
    }
    //-------------------------------------------------------------------------
    m_end;
}

void SetGaussian::ApplyOperatorF(const qid_t* qid, GridBlock* block, Field* fid) {
    //-------------------------------------------------------------------------
    real_t        pos[3];
    real_t        oo_eps2 = 1.0 / (sigma_ * sigma_);
    const real_t* xyz     = block->xyz();
    const real_t* hgrid   = block->hgrid();
   
    for (sid_t ida = 0; ida < fid->lda(); ida++) {
        real_p data = block->data(fid, ida);

        for (int i2 = 0; i2 < M_N; i2++) {
            for (int i1 = 0; i1 < M_N; i1++) {
                for (int i0 = 0; i0 < M_N; i0++) {
                    // get the position
                    m_pos(pos, i0, i1, i2, hgrid, xyz);
                    // compute the gaussian
                    const real_t x0   = pos[0] - center_[0];
                    const real_t x1   = pos[1] - center_[1];
                    const real_t x2   = pos[2] - center_[2];
                    const real_t rho2 = x0 * x0 + x1 * x1 + x2 * x2;
                    
                    data[m_idx(i0, i1, i2)] = oo_eps2 * std::exp(-rho2*oo_eps2);
                }
            }
        }
    }
    //-------------------------------------------------------------------------
}

SetAbs::SetAbs(real_t alpha[3], real_t center[3]) {
    m_begin;
    //-------------------------------------------------------------------------
    for (int id = 0; id < 3; id++) {
        center_[id] = center[id];
        alpha_[id]  = alpha[id];
    }
    //-------------------------------------------------------------------------
    m_end;
}

void SetAbs::ApplyOperatorF(const qid_t* qid, GridBlock* block, Field* fid) {
    //-------------------------------------------------------------------------
    real_t        pos[3];
    const real_t* xyz     = block->xyz();
    const real_t* hgrid   = block->hgrid();
   
    for (sid_t ida = 0; ida < fid->lda(); ida++) {
        real_p data = block->data(fid, ida);

        for (int i2 = 0; i2 < M_N; i2++) {
            for (int i1 = 0; i1 < M_N; i1++) {
                for (int i0 = 0; i0 < M_N; i0++) {
                    // get the position
                    m_pos(pos, i0, i1, i2, hgrid, xyz);
                    // compute the gaussian
                    const real_t x0   = pos[0] - center_[0];
                    const real_t x1   = pos[1] - center_[1];
                    const real_t x2   = pos[2] - center_[2];
                    
                    data[m_idx(i0, i1, i2)] = alpha_[0] * std::fabs(x0) +  alpha_[1] * std::fabs(x1) +  alpha_[2] * std::fabs(x2);
                }
            }
        }
    }
    //-------------------------------------------------------------------------
}



SetJump::SetJump(real_t alpha[3], real_t center[3]) {
    m_begin;
    //-------------------------------------------------------------------------
    for (int id = 0; id < 3; id++) {
        center_[id] = center[id];
        alpha_[id]  = alpha[id];
    }
    //-------------------------------------------------------------------------
    m_end;
}

void SetJump::ApplyOperatorF(const qid_t* qid, GridBlock* block, Field* fid) {
    //-------------------------------------------------------------------------
    real_t        pos[3];
    const real_t* xyz   = block->xyz();
    const real_t* hgrid = block->hgrid();

    for (sid_t ida = 0; ida < fid->lda(); ida++) {
        real_p data = block->data(fid, ida);

        for (int i2 = 0; i2 < M_N; i2++) {
            for (int i1 = 0; i1 < M_N; i1++) {
                for (int i0 = 0; i0 < M_N; i0++) {
                    // get the position
                    m_pos(pos, i0, i1, i2, hgrid, xyz);
                    // compute the gaussian
                    const real_t x0 = pos[0] - center_[0];
                    const real_t x1 = pos[1] - center_[1];
                    const real_t x2 = pos[2] - center_[2];

                    data[m_idx(i0, i1, i2)] = alpha_[0] * (x0 >= 0.0) + alpha_[1] * (x1 >= 0.0) + alpha_[2] * (x2 >= 0.0);
                }
            }
        }
    }
    //-------------------------------------------------------------------------
}