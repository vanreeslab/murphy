#include "setvalues.hpp"

#include <cmath>


SetGaussian::SetGaussian(real_t sigma, real_t center[3]) {
    m_begin;
    //-------------------------------------------------------------------------
    sigma_ = sigma;
    for (int id = 0; id < 3; id++) {
        center_[id] = center[id];
    }
    //-------------------------------------------------------------------------
    m_end;
}

void SetGaussian::ApplyOpF(const qid_t* qid, GridBlock* block, Field* fid) {
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

                    data[m_idx(i0, i1, i2)] = oo_eps2 * std::exp(-rho2 * oo_eps2);
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

void SetAbs::ApplyOpF(const qid_t* qid, GridBlock* block, Field* fid) {
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

                    data[m_idx(i0, i1, i2)] = alpha_[0] * std::fabs(x0) + alpha_[1] * std::fabs(x1) + alpha_[2] * std::fabs(x2);
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

void SetJump::ApplyOpF(const qid_t* qid, GridBlock* block, Field* fid) {
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

SetSinus::SetSinus(real_t length[3], real_t freq[3]) {
    m_begin;
    //-------------------------------------------------------------------------
    for (int id = 0; id < 3; id++) {
        length_[id] = length[id];
        freq_[id]   = freq[id];
    }
    //-------------------------------------------------------------------------
    m_end;
}

void SetSinus::ApplyOpF(const qid_t* qid, GridBlock* block, Field* fid) {
    //-------------------------------------------------------------------------
    real_t        pos[3];
    const real_t* xyz     = block->xyz();
    const real_t* hgrid   = block->hgrid();
    const real_t  vol[3]  = {hgrid[1] * hgrid[2], hgrid[0] * hgrid[2], hgrid[0] * hgrid[1]};
    const real_t  fact[3] = {2.0 * M_PI * freq_[0] / length_[0], 2.0 * M_PI * freq_[1] / length_[1], 2.0 * M_PI * freq_[2] / length_[2]};

    for (sid_t ida = 0; ida < fid->lda(); ida++) {
        real_p data = block->data(fid, ida);

        for (int i2 = 0; i2 < M_N; i2++) {
            for (int i1 = 0; i1 < M_N; i1++) {
                for (int i0 = 0; i0 < M_N; i0++) {
                    // get the position
                    m_pos(pos, i0, i1, i2, hgrid, xyz);

                    data[m_idx(i0, i1, i2)] = (-cos((pos[0] + hgrid[0] * 0.5) * fact[0]) + cos((pos[0] - hgrid[0] * 0.5) * fact[0])) / (hgrid[0] * fact[0]) *
                                              (-cos((pos[1] + hgrid[1] * 0.5) * fact[1]) + cos((pos[1] - hgrid[1] * 0.5) * fact[1])) / (hgrid[1] * fact[1]) *
                                              (-cos((pos[2] + hgrid[2] * 0.5) * fact[2]) + cos((pos[2] - hgrid[2] * 0.5) * fact[2])) / (hgrid[2] * fact[2]);
                }
            }
        }
    }
    //-------------------------------------------------------------------------
}

SetCosinus::SetCosinus(real_t length[3], real_t freq[3]) {
    m_begin;
    //-------------------------------------------------------------------------
    for (int id = 0; id < 3; id++) {
        length_[id] = length[id];
        freq_[id]   = freq[id];
    }
    //-------------------------------------------------------------------------
    m_end;
}

void SetCosinus::ApplyOpF(const qid_t* qid, GridBlock* block, Field* fid) {
    //-------------------------------------------------------------------------
    real_t        pos[3];
    const real_t* xyz     = block->xyz();
    const real_t* hgrid   = block->hgrid();
    const real_t  vol[3]  = {hgrid[1] * hgrid[2], hgrid[0] * hgrid[2], hgrid[0] * hgrid[1]};
    const real_t  fact[3] = {2.0 * M_PI * freq_[0] / length_[0], 2.0 * M_PI * freq_[1] / length_[1], 2.0 * M_PI * freq_[2] / length_[2]};

    for (sid_t ida = 0; ida < fid->lda(); ida++) {
        real_p data = block->data(fid, ida);

        for (int i2 = 0; i2 < M_N; i2++) {
            for (int i1 = 0; i1 < M_N; i1++) {
                for (int i0 = 0; i0 < M_N; i0++) {
                    // get the position
                    m_pos(pos, i0, i1, i2, hgrid, xyz);

                    data[m_idx(i0, i1, i2)] = (sin((pos[0] + hgrid[0] * 0.5) * fact[0]) - sin((pos[0] - hgrid[0] * 0.5) * fact[0])) / hgrid[0] / fact[0] +
                                              (sin((pos[1] + hgrid[1] * 0.5) * fact[1]) - sin((pos[1] - hgrid[1] * 0.5) * fact[1])) / hgrid[1] / fact[1] +
                                              (sin((pos[2] + hgrid[2] * 0.5) * fact[2]) - sin((pos[2] - hgrid[2] * 0.5) * fact[2])) / hgrid[2] / fact[2];
                    // data[m_idx(i0,i1,i2)] = 0.25*(pow(pos[0]+hgrid[0]*0.5,4)-pow(pos[0]-hgrid[0]*0.5,4)) / (hgrid[0]);
                    // data[m_idx(i0,i1,i2)] = 1.0/6.0*(pow(pos[0]+hgrid[0]*0.5,6)-pow(pos[0]-hgrid[0]*0.5,6)) / (hgrid[0]);
                }
            }
        }
    }
    //-------------------------------------------------------------------------
}

SetPolynom::SetPolynom(lid_t degree[3], real_t direction[3]) {
    m_begin;
    //-------------------------------------------------------------------------
    for (int id = 0; id < 3; id++) {
        deg_[id] = degree[id];
        dir_[id] = direction[id];
    }
    //-------------------------------------------------------------------------
    m_end;
}
void SetPolynom::ApplyOpF(const qid_t* qid, GridBlock* block, Field* fid) {
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

                    data[m_idx(i0, i1, i2)] = dir_[0] / (deg_[0] + 1.0) * (pow(pos[0] + hgrid[0] * 0.5, deg_[0] + 1.0) - pow(pos[0] - hgrid[0] * 0.5, deg_[0] + 1.0)) / hgrid[0] +
                                              dir_[1] / (deg_[1] + 1.0) * (pow(pos[1] + hgrid[1] * 0.5, deg_[1] + 1.0) - pow(pos[1] - hgrid[1] * 0.5, deg_[1] + 1.0)) / hgrid[1] +
                                              dir_[2] / (deg_[2] + 1.0) * (pow(pos[2] + hgrid[2] * 0.5, deg_[2] + 1.0) - pow(pos[2] - hgrid[2] * 0.5, deg_[2] + 1.0)) / hgrid[2];
                }
            }
        }
    }
    //-------------------------------------------------------------------------
}

SetExponential::SetExponential(real_t center[3], real_t sigma[3], real_t alpha) {
    m_begin;
    //-------------------------------------------------------------------------
    for (int id = 0; id < 3; id++) {
        center_[id] = center[id];
        sigma_[id]  = sigma[id];
    }
    alpha_ = alpha;
    //-------------------------------------------------------------------------
    m_end;
}

void SetExponential::ApplyOpF(const qid_t* qid, GridBlock* block, Field* fid) {
    //-------------------------------------------------------------------------
    real_t        pos[3];
    const real_t* xyz   = block->xyz();
    const real_t* hgrid = block->hgrid();

    real_t sigma     = sqrt(sigma_[0] * sigma_[0] + sigma_[1] * sigma_[1] + sigma_[2] * sigma_[2]);
    real_t oo_sigma2 = (sigma > 0.0) ? 1.0 / (sigma * sigma) : 0.0;
    real_t fact = alpha_* sqrt(2.0 / M_PI) / (4.0 * M_PI * pow(sigma, 3.0)); // see Wincky encyclopedia

    for (sid_t ida = 0; ida < fid->lda(); ida++) {
        real_p data = block->data(fid, ida);

        for (int i2 = 0; i2 < M_N; i2++) {
            for (int i1 = 0; i1 < M_N; i1++) {
                for (int i0 = 0; i0 < M_N; i0++) {
                    // get the position
                    m_pos(pos, i0, i1, i2, hgrid, xyz);
                    // compute the gaussian
                    const real_t x0         = (sigma_[0] > 0) ? pos[0] - center_[0] : 0.0;
                    const real_t x1         = (sigma_[1] > 0) ? pos[1] - center_[1] : 0.0;
                    const real_t x2         = (sigma_[2] > 0) ? pos[2] - center_[2] : 0.0;
                    const real_t rho2       = (x0 * x0 + x1 * x1 + x2 * x2) * oo_sigma2;
                    data[m_idx(i0, i1, i2)] = fact * std::exp(-rho2 * 0.5);
                }
            }
        }
    }
    //-------------------------------------------------------------------------
}

SetErf::SetErf(real_t center[3], real_t sigma[3], real_t alpha) {
    m_begin;
    //-------------------------------------------------------------------------
    for (int id = 0; id < 3; id++) {
        center_[id] = center[id];
        sigma_[id]  = sigma[id];
    }
    alpha_ = alpha;
    //-------------------------------------------------------------------------
    m_end;
}

void SetErf::ApplyOpF(const qid_t* qid, GridBlock* block, Field* fid) {
    //-------------------------------------------------------------------------
    real_t        pos[3];
    const real_t* xyz   = block->xyz();
    const real_t* hgrid = block->hgrid();

    real_t sigma     = sqrt(sigma_[0] * sigma_[0] + sigma_[1] * sigma_[1] + sigma_[2] * sigma_[2]);
    real_t oo_sigma2 = (sigma > 0.0) ? 1.0 / (sigma * sigma) : 0.0;
    real_t oo_sqrt2  = 1.0 / M_SQRT2;
    real_t fact      = alpha_ / (4.0*M_PI*sigma);  // see Wincky encyclopedia

    for (sid_t ida = 0; ida < fid->lda(); ida++) {
        real_p data = block->data(fid, ida);

        for (int i2 = 0; i2 < M_N; i2++) {
            for (int i1 = 0; i1 < M_N; i1++) {
                for (int i0 = 0; i0 < M_N; i0++) {
                    // get the position
                    m_pos(pos, i0, i1, i2, hgrid, xyz);
                    // compute the gaussian
                    const real_t x0         = (sigma_[0] > 0) ? pos[0] - center_[0] : 0.0;
                    const real_t x1         = (sigma_[1] > 0) ? pos[1] - center_[1] : 0.0;
                    const real_t x2         = (sigma_[2] > 0) ? pos[2] - center_[2] : 0.0;
                    const real_t rho2       = (x0 * x0 + x1 * x1 + x2 * x2) * oo_sigma2;
                    const real_t rho        = sqrt(rho2);
                    data[m_idx(i0, i1, i2)] = fact / rho * std::erf(rho * oo_sqrt2);
                }
            }
        }
    }
    //-------------------------------------------------------------------------
}

// SetExpoCosinus::SetExpoCosinus(real_t center[3], real_t sigma[3],real_t length[3], real_t freq[3]) {
//     m_begin;
//     //-------------------------------------------------------------------------
//     for (int id = 0; id < 3; id++) {
//         length_[id] = length[id];
//         freq_[id]   = freq[id];
//         center_[id] = center[id];
//         sigma_[id]  = sigma[id];
//     }
//     //-------------------------------------------------------------------------
//     m_end;
// }

// void SetExpoCosinus::ApplyOpF(const qid_t* qid, GridBlock* block, Field* fid) {
//     //-------------------------------------------------------------------------
//     real_t        pos[3];
//     const real_t* xyz   = block->xyz();
//     const real_t* hgrid = block->hgrid();

//     real_t sigma          = sqrt(sigma_[0] * sigma_[0] + sigma_[1] * sigma_[1] + sigma_[2] * sigma_[2]);
//     real_t oo_eps2 = (sigma > 0.0) ? 1.0 / (sigma * sigma) : 0.0;

//     for (sid_t ida = 0; ida < fid->lda(); ida++) {
//         real_p data = block->data(fid, ida);

//         for (int i2 = 0; i2 < M_N; i2++) {
//             for (int i1 = 0; i1 < M_N; i1++) {
//                 for (int i0 = 0; i0 < M_N; i0++) {
//                     // get the position
//                     m_pos(pos, i0, i1, i2, hgrid, xyz);

//                     // compute the gaussian
//                     const real_t x0   = (sigma_[0] > 0) ? pos[0] - center_[0] : 0.0;
//                     const real_t x1   = (sigma_[1] > 0) ? pos[1] - center_[1] : 0.0;
//                     const real_t x2   = (sigma_[2] > 0) ? pos[2] - center_[2] : 0.0;
//                     const real_t rho2 = x0 * x0 + x1 * x1 + x2 * x2;

//                     data[m_idx(i0, i1, i2)] = std::cos(2.0 * M_PI * pos[0] / length_[0] * freq_[0]) *
//                                               std::cos(2.0 * M_PI * pos[1] / length_[1] * freq_[1]) *
//                                               std::cos(2.0 * M_PI * pos[2] / length_[2] * freq_[2]) *
//                                               oo_eps2 * std::exp(-rho2 * oo_eps2);
//                 }
//             }
//         }
//     }
//     //-------------------------------------------------------------------------
// }