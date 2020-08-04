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
    const real_t* xyz   = block->xyz();
    const real_t* hgrid   = block->hgrid();
    const real_t  fact[3] = {2.0 * M_PI * freq_[0] / length_[0], 2.0 * M_PI * freq_[1] / length_[1], 2.0 * M_PI * freq_[2] / length_[2]};

    for (sid_t ida = 0; ida < fid->lda(); ida++) {
        real_p data = block->data(fid, ida);

        for (int i2 = 0; i2 < M_N; i2++) {
            for (int i1 = 0; i1 < M_N; i1++) {
                for (int i0 = 0; i0 < M_N; i0++) {
                    // get the position
                    m_pos(pos, i0, i1, i2, hgrid, xyz);

                    data[m_idx(i0, i1, i2)] = sin(pos[0] * fact[0]) * sin(pos[1] * fact[1]) * sin(pos[2] * fact[2]);
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
    const real_t  vol     = 1.0 / (hgrid[0] * hgrid[1] * hgrid[2]);
    const real_t  fact[3] = {2.0 * M_PI * freq_[0] / length_[0], 2.0 * M_PI * freq_[1] / length_[1], 2.0 * M_PI * freq_[2] / length_[2]};

    for (sid_t ida = 0; ida < fid->lda(); ida++) {
        real_p data = block->data(fid, ida);

        for (int i2 = 0; i2 < M_N; i2++) {
            for (int i1 = 0; i1 < M_N; i1++) {
                for (int i0 = 0; i0 < M_N; i0++) {
                    // get the position
                    m_pos(pos, i0, i1, i2, hgrid, xyz);
                    // data[m_idx(i0, i1, i2)] = vol * ((sin((pos[0] + hgrid[0] * 0.5) * fact[0]) - sin((pos[0] - hgrid[0] * 0.5) * fact[0])) / fact[0] *
                    //                                  (sin((pos[1] + hgrid[1] * 0.5) * fact[1]) - sin((pos[1] - hgrid[1] * 0.5) * fact[1])) / fact[1] *
                    //                                  (sin((pos[2] + hgrid[2] * 0.5) * fact[2]) - sin((pos[2] - hgrid[2] * 0.5) * fact[2])) / fact[2]);
                    data[m_idx(i0, i1, i2)] = cos(pos[0] * fact[0]) *
                                              cos(pos[1] * fact[1]) *
                                              cos(pos[2] * fact[2]);
                }
            }
        }
    }
    //-------------------------------------------------------------------------
}

SetLaplaCosinus::SetLaplaCosinus(real_t length[3], real_t freq[3]) {
    m_begin;
    //-------------------------------------------------------------------------
    for (int id = 0; id < 3; id++) {
        length_[id] = length[id];
        freq_[id]   = freq[id];
    }
    //-------------------------------------------------------------------------
    m_end;
}

void SetLaplaCosinus::ApplyOpF(const qid_t* qid, GridBlock* block, Field* fid) {
    //-------------------------------------------------------------------------
    real_t        pos[3];
    const real_t* xyz     = block->xyz();
    const real_t* hgrid   = block->hgrid();
    const real_t  fact[3] = {2.0 * M_PI * freq_[0] / length_[0], 2.0 * M_PI * freq_[1] / length_[1], 2.0 * M_PI * freq_[2] / length_[2]};
    const real_t  vol     = -(fact[0] * fact[0] + fact[1] * fact[1] + fact[2] * fact[2]) / (hgrid[0] * hgrid[1] * hgrid[2]);

    for (sid_t ida = 0; ida < fid->lda(); ida++) {
        real_p data = block->data(fid, ida);
        for (int i2 = 0; i2 < M_N; i2++) {
            for (int i1 = 0; i1 < M_N; i1++) {
                for (int i0 = 0; i0 < M_N; i0++) {
                    // get the position
                    m_pos(pos, i0, i1, i2, hgrid, xyz);
                    data[m_idx(i0, i1, i2)] = vol * ((sin((pos[0] + hgrid[0] * 0.5) * fact[0]) - sin((pos[0] - hgrid[0] * 0.5) * fact[0])) / fact[0] *
                                                     (sin((pos[1] + hgrid[1] * 0.5) * fact[1]) - sin((pos[1] - hgrid[1] * 0.5) * fact[1])) / fact[1] *
                                                     (sin((pos[2] + hgrid[2] * 0.5) * fact[2]) - sin((pos[2] - hgrid[2] * 0.5) * fact[2])) / fact[2]);
                }
            }
        }
    }
    //-------------------------------------------------------------------------
}

SetPolynom::SetPolynom(const lid_t degree[3], const real_t direction[3]) : SetPolynom(degree, direction, false) {
    m_begin;
    //-------------------------------------------------------------------------
    // do nothing
    //-------------------------------------------------------------------------
    m_end;
}

SetPolynom::SetPolynom(const lid_t degree[3], const real_t direction[3], const bool extend) {
    m_begin;
    //-------------------------------------------------------------------------
    extend_ = extend;
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

    const lid_t start = (extend_) ? (-M_GS) : (0);
    const lid_t end   = (extend_) ? (M_N + M_GS) : (M_N);

    for (sid_t ida = 0; ida < fid->lda(); ida++) {
        real_p data = block->data(fid, ida);

        for (int i2 = start; i2 < end; i2++) {
            for (int i1 = start; i1 < end; i1++) {
                for (int i0 = start; i0 < end; i0++) {
                    // get the position
                    m_pos(pos, i0, i1, i2, hgrid, xyz);

                    data[m_idx(i0, i1, i2)] = dir_[0] * pow(pos[0], deg_[0]) + dir_[1] * pow(pos[1], deg_[1]) + dir_[2] * pow(pos[2], deg_[2]);
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
    for (lda_t id = 0; id < 3; id++) {
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

SetVortexRing::SetVortexRing(const lda_t normal, const real_t center[3], const real_t sigma, const real_t radius) {
    m_begin;
    //-------------------------------------------------------------------------
    normal_ = normal;
    sigma_  = sigma;
    radius_ = radius;
    for (lda_t id = 0; id < 3; ++id) {
        center_[id] = center[id];
    }
    //-------------------------------------------------------------------------
    m_end;
}

void SetVortexRing::ApplyOpF(const qid_t* qid, GridBlock* block, Field* fid) {
    //-------------------------------------------------------------------------
    real_t        pos[3];
    const real_t* xyz   = block->xyz();
    const real_t* hgrid = block->hgrid();

    const real_t oo_sigma2   = 1.0 / (sigma_ * sigma_);
    const real_t oo_pisigma2 = 1.0 / (M_PI * sigma_ * sigma_);

    // compute the normal direction as being the z one and the two other as x and y
    const lda_t idx = (normal_ + 1) % 3;
    const lda_t idy = (normal_ + 2) % 3;
    const lda_t idz = normal_;
    // get the pointers correct
    data_ptr wx = block->data(fid, idx);
    data_ptr wy = block->data(fid, idy);
    data_ptr wz = block->data(fid, idz);

    for (int i2 = (-M_GS); i2 < (M_N + M_GS); i2++) {
        for (int i1 = (-M_GS); i1 < (M_N + M_GS); i1++) {
            for (int i0 = (-M_GS); i0 < (M_N + M_GS); i0++) {
                // get the position
                m_pos(pos, i0, i1, i2, hgrid, xyz);
                // wrt to the center
                const real_t alpha = atan2(pos[idy] - center_[idy], pos[idx] - center_[idx]);
                const real_t x     = pos[idx] - (center_[idx] + radius_ * cos(alpha));
                const real_t y     = pos[idy] - (center_[idy] + radius_ * sin(alpha));
                const real_t z     = pos[idz] - (center_[idz]);
                // get the local coords

                // const real_t r_plane = sqrt(x * x + y * y) - radius_;
                const real_t r2 = (x * x + y * y + z * z);
                const real_t vort = oo_pisigma2 * exp(-r2 * oo_sigma2);
                // const real_t vort = oo_pisigma2 * exp(1.0 - 1.0 / (1.0 - r2 * oo_sigma2));

                wx[m_idx(i0, i1, i2)] = -vort * sin(alpha);  //  -vort * cos(alpha);
                wy[m_idx(i0, i1, i2)] = vort * cos(alpha);   // +vort * sin(alpha);
                wz[m_idx(i0, i1, i2)] = 0.0;
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
