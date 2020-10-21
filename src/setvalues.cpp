#include "setvalues.hpp"

#include <cmath>

/**
 * @brief Set the Value and defines the range which is assigned to the value, including ghost points or not
 * 
 * @param interp the Interpolator object used to retrieve the objects
 * @param profiler the profiler used to time the operations
 */
SetValue::SetValue(const InterpolatingWavelet* interp) : BlockOperator(interp) {
    m_begin;
    //-------------------------------------------------------------------------
    // m_profCreate(profiler, "SetValue");
    //-------------------------------------------------------------------------
    m_end;
}

void SetValue::operator()(const ForestGrid* grid, Field* field) {
    m_begin;
    //-------------------------------------------------------------------------
    // get the span of ida
    ida_start_ = 0;
    ida_end_   = field->lda();
    // go for it
    DoOpTree(this, &SetValue::FillGridBlock, grid, field);
    // update the ghost status
    m_verb("setting the ghosts of %s to false", field->name().c_str());
    field->ghost_status(do_ghost_);
    //-------------------------------------------------------------------------
    m_end;
}

void SetValue::operator()(const ForestGrid* grid, Field* field, const lda_t ida) {
    m_begin;
    //-------------------------------------------------------------------------
    // get the span of ida
    ida_start_ = ida;
    ida_end_   = ida + 1;
    // go for it
    DoOpTree(this, &SetValue::FillGridBlock, grid, field);
    // update the ghost status
    m_verb("setting the ghosts of %s to false", field->name().c_str());
    field->ghost_status(do_ghost_);
    //-------------------------------------------------------------------------
    m_end;
}

//=====================================================================================================
SetAbs::SetAbs(const real_t alpha[3], const real_t center[3]) : SetAbs(alpha, center, nullptr) {}
SetAbs::SetAbs(const real_t alpha[3], const real_t center[3], const InterpolatingWavelet* interp) : SetValue(interp) {
    m_begin;
    //-------------------------------------------------------------------------
    for (lda_t id = 0; id < 3; id++) {
        center_[id] = center[id];
        alpha_[id]  = alpha[id];
    }
    //-------------------------------------------------------------------------
    m_end;
}

void SetAbs::FillGridBlock(const qid_t* qid, GridBlock* block, Field* fid) {
    //-------------------------------------------------------------------------
    real_t        pos[3];
    const real_t* xyz   = block->xyz();
    const real_t* hgrid = block->hgrid();

    for (lda_t ida = ida_start_; ida < ida_end_; ida++) {
        real_p data = block->data(fid, ida);
        for (lid_t i2 = start_; i2 < end_; i2++) {
            for (lid_t i1 = start_; i1 < end_; i1++) {
                for (lid_t i0 = start_; i0 < end_; i0++) {
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

//=====================================================================================================
SetSinus::SetSinus(const real_t length[3], const real_t freq[3], const real_t alpha[3]) : SetSinus(length, freq, alpha, nullptr) {}
SetSinus::SetSinus(const real_t length[3], const real_t freq[3], const real_t alpha[3], const InterpolatingWavelet* interp) : SetValue(interp) {
    m_begin;
    //-------------------------------------------------------------------------
    for (lda_t id = 0; id < 3; id++) {
        length_[id] = length[id];
        freq_[id]   = freq[id];
        alpha_[id]  = alpha[id];
    }
    //-------------------------------------------------------------------------
    m_end;
}

void SetSinus::FillGridBlock(const qid_t* qid, GridBlock* block, Field* fid) {
    //-------------------------------------------------------------------------
    real_t        pos[3];
    const real_t* xyz     = block->xyz();
    const real_t* hgrid   = block->hgrid();
    const real_t  fact[3] = {2.0 * M_PI * freq_[0] / length_[0], 2.0 * M_PI * freq_[1] / length_[1], 2.0 * M_PI * freq_[2] / length_[2]};

    for (lda_t ida = ida_start_; ida < ida_end_; ida++) {
        real_p data = block->data(fid, ida);
        for (lid_t i2 = start_; i2 < end_; i2++) {
            for (lid_t i1 = start_; i1 < end_; i1++) {
                for (lid_t i0 = start_; i0 < end_; i0++) {
                    // get the position
                    m_pos(pos, i0, i1, i2, hgrid, xyz);
                    data[m_idx(i0, i1, i2)] = alpha_[0] * sin(pos[0] * fact[0]) +
                                              alpha_[1] * sin(pos[1] * fact[1]) +
                                              alpha_[2] * sin(pos[2] * fact[2]);
                }
            }
        }
    }
    //-------------------------------------------------------------------------
}

//=====================================================================================================
SetCosinus::SetCosinus(const real_t length[3], const real_t freq[3], const real_t alpha[3]) : SetCosinus(length, freq, alpha, nullptr) {}
SetCosinus::SetCosinus(const real_t length[3], const real_t freq[3], const real_t alpha[3], const InterpolatingWavelet* interp) : SetValue(interp) {
    m_begin;
    //-------------------------------------------------------------------------
    for (int id = 0; id < 3; id++) {
        length_[id] = length[id];
        freq_[id]   = freq[id];
        alpha_[id]  = alpha[id];
    }
    //-------------------------------------------------------------------------
    m_end;
}

void SetCosinus::FillGridBlock(const qid_t* qid, GridBlock* block, Field* fid) {
    //-------------------------------------------------------------------------
    real_t        pos[3];
    const real_t* xyz     = block->xyz();
    const real_t* hgrid   = block->hgrid();
    const real_t  vol     = 1.0 / (hgrid[0] * hgrid[1] * hgrid[2]);
    const real_t  fact[3] = {2.0 * M_PI * freq_[0] / length_[0], 2.0 * M_PI * freq_[1] / length_[1], 2.0 * M_PI * freq_[2] / length_[2]};

    for (lda_t ida = ida_start_; ida < ida_end_; ida++) {
        real_p data = block->data(fid, ida);
        for (lid_t i2 = start_; i2 < end_; i2++) {
            for (lid_t i1 = start_; i1 < end_; i1++) {
                for (lid_t i0 = start_; i0 < end_; i0++) {
                    // get the position
                    m_pos(pos, i0, i1, i2, hgrid, xyz);
                    data[m_idx(i0, i1, i2)] = alpha_[0] * cos(pos[0] * fact[0]) +
                                              alpha_[1] * cos(pos[1] * fact[1]) +
                                              alpha_[2] * cos(pos[2] * fact[2]);
                }
            }
        }
    }
    //-------------------------------------------------------------------------
}

//=====================================================================================================
SetPolynom::SetPolynom(const lid_t degree[3], const real_t direction[3], const real_t shift[3]) : SetPolynom(degree, direction, shift, nullptr) {}
SetPolynom::SetPolynom(const lid_t degree[3], const real_t direction[3], const real_t shift[3], const InterpolatingWavelet* interp) : SetValue(interp) {
    m_begin;
    //-------------------------------------------------------------------------
    for (lda_t id = 0; id < 3; id++) {
        deg_[id]   = degree[id];
        dir_[id]   = direction[id];
        shift_[id] = shift[id];
    }
    //-------------------------------------------------------------------------
    m_end;
}
void SetPolynom::FillGridBlock(const qid_t* qid, GridBlock* block, Field* fid) {
    //-------------------------------------------------------------------------
    real_t        pos[3];
    const real_t* xyz   = block->xyz();
    const real_t* hgrid = block->hgrid();

    for (lda_t ida = ida_start_; ida < ida_end_; ida++) {
        real_p data = block->data(fid, ida);
        for (lid_t i2 = start_; i2 < end_; i2++) {
            for (lid_t i1 = start_; i1 < end_; i1++) {
                for (lid_t i0 = start_; i0 < end_; i0++) {
                    // get the position
                    m_pos(pos, i0, i1, i2, hgrid, xyz);

                    data[m_idx(i0, i1, i2)] = dir_[0] * pow(pos[0] - shift_[0], deg_[0]) +
                                              dir_[1] * pow(pos[1] - shift_[1], deg_[1]) +
                                              dir_[2] * pow(pos[2] - shift_[2], deg_[2]);
                }
            }
        }
    }
    //-------------------------------------------------------------------------
}

//=====================================================================================================
SetExponential::SetExponential(const real_t center[3], const real_t sigma[3], const real_t alpha) : SetExponential(center, sigma, alpha, nullptr) {}
SetExponential::SetExponential(const real_t center[3], const real_t sigma[3], const real_t alpha, const InterpolatingWavelet* interp) : SetValue(interp) {
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

void SetExponential::FillGridBlock(const qid_t* qid, GridBlock* block, Field* fid) {
    //-------------------------------------------------------------------------
    real_t        pos[3];
    const real_t* xyz   = block->xyz();
    const real_t* hgrid = block->hgrid();

    real_t sigma     = sqrt(sigma_[0] * sigma_[0] + sigma_[1] * sigma_[1] + sigma_[2] * sigma_[2]);
    real_t oo_sigma2 = (sigma > 0.0) ? 1.0 / (sigma * sigma) : 0.0;
    real_t fact      = alpha_ * sqrt(1.0 / M_PI * oo_sigma2);

    for (lda_t ida = ida_start_; ida < ida_end_; ida++) {
        real_p data = block->data(fid, ida);
        for (lid_t i2 = start_; i2 < end_; i2++) {
            for (lid_t i1 = start_; i1 < end_; i1++) {
                for (lid_t i0 = start_; i0 < end_; i0++) {
                    // get the position
                    m_pos(pos, i0, i1, i2, hgrid, xyz);
                    // compute the gaussian
                    const real_t rhox = (sigma_[0] > 0) ? ((pos[0] - center_[0]) / sigma_[0]) : 0.0;
                    const real_t rhoy = (sigma_[1] > 0) ? ((pos[1] - center_[1]) / sigma_[1]) : 0.0;
                    const real_t rhoz = (sigma_[2] > 0) ? ((pos[2] - center_[2]) / sigma_[2]) : 0.0;
                    const real_t rho  = rhox * rhox + rhoy * rhoy + rhoz * rhoz;

                    data[m_idx(i0, i1, i2)] = fact * std::exp(-rho);
                }
            }
        }
    }
    //-------------------------------------------------------------------------
}

//=====================================================================================================
SetErf::SetErf(const real_t center[3], const real_t sigma[3], const real_t alpha) : SetErf(center, sigma, alpha, nullptr) {}
SetErf::SetErf(const real_t center[3], const real_t sigma[3], const real_t alpha, const InterpolatingWavelet* interp) : SetValue(interp) {
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

void SetErf::FillGridBlock(const qid_t* qid, GridBlock* block, Field* fid) {
    //-------------------------------------------------------------------------
    real_t        pos[3];
    const real_t* xyz   = block->xyz();
    const real_t* hgrid = block->hgrid();

    real_t sigma     = sqrt(sigma_[0] * sigma_[0] + sigma_[1] * sigma_[1] + sigma_[2] * sigma_[2]);
    real_t oo_sigma2 = (sigma > 0.0) ? 1.0 / (sigma * sigma) : 0.0;
    real_t oo_sqrt2  = 1.0 / M_SQRT2;
    real_t fact      = alpha_ / (4.0 * M_PI * sigma);  // see Wincky encyclopedia

    for (lda_t ida = ida_start_; ida < ida_end_; ida++) {
        real_p data = block->data(fid, ida);
        for (lid_t i2 = start_; i2 < end_; i2++) {
            for (lid_t i1 = start_; i1 < end_; i1++) {
                for (lid_t i0 = start_; i0 < end_; i0++) {
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

//=====================================================================================================
SetVortexRing::SetVortexRing(const lda_t normal, const real_t center[3], const real_t sigma, const real_t radius) : SetVortexRing(normal, center, sigma, radius, nullptr) {}
SetVortexRing::SetVortexRing(const lda_t normal, const real_t center[3], const real_t sigma, const real_t radius, const InterpolatingWavelet* interp) : SetValue(interp) {
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

void SetVortexRing::FillGridBlock(const qid_t* qid, GridBlock* block, Field* fid) {
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

    for (lid_t i2 = start_; i2 < end_; i2++) {
        for (lid_t i1 = start_; i1 < end_; i1++) {
            for (lid_t i0 = start_; i0 < end_; i0++) {
                // get the position
                m_pos(pos, i0, i1, i2, hgrid, xyz);
                // wrt to the center
                const real_t alpha = atan2(pos[idy] - center_[idy], pos[idx] - center_[idx]);
                const real_t x     = pos[idx] - (center_[idx] + radius_ * cos(alpha));
                const real_t y     = pos[idy] - (center_[idy] + radius_ * sin(alpha));
                const real_t z     = pos[idz] - (center_[idz]);
                // get the local coords

                // const real_t r_plane = sqrt(x * x + y * y) - radius_;
                const real_t r2   = (x * x + y * y + z * z);
                const real_t rho2 = r2 * oo_sigma2;
                const real_t vort = oo_pisigma2 * exp(-rho2);
                // const real_t vort = oo_pisigma2 * exp(1.0 - 1.0 / (1.0 - r2 * oo_sigma2));
                wx[m_idx(i0, i1, i2)] = -vort * sin(alpha);  //  -vort * cos(alpha);
                wy[m_idx(i0, i1, i2)] = vort * cos(alpha);   // +vort * sin(alpha);
                wz[m_idx(i0, i1, i2)] = 0.0;
            }
        }
    }
    //-------------------------------------------------------------------------
}

//=====================================================================================================
SetCompactVortexRing::SetCompactVortexRing(const lda_t normal, const real_t center[3], const real_t sigma, const real_t radius, const real_t cutoff) : SetCompactVortexRing(normal, center, sigma, radius, cutoff, nullptr) {}
SetCompactVortexRing::SetCompactVortexRing(const lda_t normal, const real_t center[3], const real_t sigma, const real_t radius, const real_t cutoff, const InterpolatingWavelet* interp) : SetValue(interp) {
    m_begin;
    //-------------------------------------------------------------------------
    normal_ = normal;
    sigma_  = sigma;
    radius_ = radius;
    cutoff_ = cutoff;
    for (lda_t id = 0; id < 3; ++id) {
        center_[id] = center[id];
    }
    //-------------------------------------------------------------------------
    m_end;
}

void SetCompactVortexRing::FillGridBlock(const qid_t* qid, GridBlock* block, Field* fid) {
    //-------------------------------------------------------------------------
    real_t        pos[3];
    const real_t* xyz   = block->xyz();
    const real_t* hgrid = block->hgrid();

    const real_t oo_sigma2   = 1.0 / (sigma_ * sigma_);
    const real_t oo_R2       = 1.0 / (cutoff_ * cutoff_);
    const real_t oo_pisigma2 = 1.0 / (M_PI * sigma_ * sigma_);

    // compute the normal direction as being the z one and the two other as x and y
    const lda_t idx = (normal_ + 1) % 3;
    const lda_t idy = (normal_ + 2) % 3;
    const lda_t idz = normal_;
    // get the pointers correct
    data_ptr wx = block->data(fid, idx);
    data_ptr wy = block->data(fid, idy);
    data_ptr wz = block->data(fid, idz);

    for (lid_t i2 = start_; i2 < end_; i2++) {
        for (lid_t i1 = start_; i1 < end_; i1++) {
            for (lid_t i0 = start_; i0 < end_; i0++) {
                // get the position
                m_pos(pos, i0, i1, i2, hgrid, xyz);
                // wrt to the center
                const real_t alpha = atan2(pos[idy] - center_[idy], pos[idx] - center_[idx]);
                const real_t x     = pos[idx] - (center_[idx] + radius_ * cos(alpha));
                const real_t y     = pos[idy] - (center_[idy] + radius_ * sin(alpha));
                const real_t z     = pos[idz] - (center_[idz]);
                // get the local coords

                // const real_t r_plane = sqrt(x * x + y * y) - radius_;
                const real_t r2    = (x * x + y * y + z * z);
                const real_t rho2  = r2 * oo_sigma2;
                const real_t rhoR2 = r2 * oo_R2;

                if (rhoR2 < 1.0) {
                    const real_t vort = oo_pisigma2 * exp(-rho2 / (1.0 - rhoR2));
                    // const real_t vort = oo_pisigma2 * exp(1.0 - 1.0 / (1.0 - r2 * oo_sigma2));
                    wx[m_idx(i0, i1, i2)] = -vort * sin(alpha);  //  -vort * cos(alpha);
                    wy[m_idx(i0, i1, i2)] = vort * cos(alpha);   // +vort * sin(alpha);
                    wz[m_idx(i0, i1, i2)] = 0.0;

                } else {
                    wx[m_idx(i0, i1, i2)] = 0.0;
                    wy[m_idx(i0, i1, i2)] = 0.0;
                    wz[m_idx(i0, i1, i2)] = 0.0;
                }
            }
        }
    }
    //-------------------------------------------------------------------------
}
