#include "setvalues.hpp"

/**
 * @brief returns the value of an exponential blob between 0 and 1 (given its sigma and its center)
 * 
 */
lambda_t<real_t, const real_t[], const real_t[], const real_t>
    scalar_exp = [](const real_t pos[3], const real_t center[3], const real_t sigma) -> real_t {
    const real_t fact = 1.0;  //sqrt(M_PI*sigma*sigma) / sqrt(M_PI * diff_sigma2 );
    // compute the gaussian
    const real_t rhox = (pos[0] - center[0]) / sigma;
    const real_t rhoy = (pos[1] - center[1]) / sigma;
    const real_t rhoz = (pos[2] - center[2]) / sigma;
    const real_t rho  = rhox * rhox + rhoy * rhoy + rhoz * rhoz;
    return fact * std::exp(-rho);
};

/**
 * @brief diffusive exponential
 * 
 * the integral of the values is sqrt(M_PI * sigma^2)
 */
lambda_t<real_t, const real_t[], const real_t[], const real_t, const real_t, const real_t>
    scalar_diff_exp = [](const real_t pos[3], const real_t center[3], const real_t sigma, const real_t time, const real_t nu) -> real_t {
    const real_t sigma2      = sigma * sigma;
    const real_t diff_sigma2 = sigma2 + 4.0 * nu * time;
    const real_t fact        = sqrt(M_PI * sigma2) / sqrt(M_PI * diff_sigma2);

    // compute the gaussian
    const real_t rhox = (pos[0] - center[0]);
    const real_t rhoy = (pos[1] - center[1]);
    const real_t rhoz = (pos[2] - center[2]);
    const real_t rho  = (rhox * rhox + rhoy * rhoy + rhoz * rhoz) / diff_sigma2;
    return fact * std::exp(-rho);
};

/**
 * @brief returns the value of an compact exponential blob between 0 and 1 (given its sigma and its center)
 * 
 */
lambda_t<real_t, const real_t[], const real_t[], const real_t, const real_t>
    scalar_compact_exp = [](const real_t pos[3], const real_t center[3], const real_t sigma, const real_t beta) -> real_t {
    const real_t gamma     = beta * sigma;
    const real_t oo_sigma2 = 1.0 / (sigma * sigma);
    const real_t oo_gamma2 = 1.0 / (gamma * gamma);
    // compute the gaussian
    const real_t x      = (pos[0] - center[0]);
    const real_t y      = (pos[1] - center[1]);
    const real_t z      = (pos[2] - center[2]);
    const real_t rad_sq = x * x + y * y + z * z;
    const real_t rho1   = rad_sq * oo_sigma2;
    const real_t rho2   = rad_sq * oo_gamma2;
    const real_t vort   = (rho2 < 1.0) ? (exp(-rho1 / (1.0 - rho2))) : 0.0;
    return vort;
};

/**
 * @brief returns the value of a scalar exponential tube aligned with the normal direction
 * 
 */
lambda_t<real_t, const real_t[], const real_t[], const real_t, const lda_t>
    scalar_tube = [](const real_t pos[3], const real_t center[3], const real_t sigma, const lda_t normal) -> real_t {
    const lda_t  idx       = (normal + 1) % 3;
    const lda_t  idy       = (normal + 2) % 3;
    const real_t oo_sigma2 = 1.0 / (sigma * sigma);
    // compute the gaussian
    const real_t rad1 = pow(pos[idx] - (center[idx]), 2) + pow(pos[idy] - center[idy], 2);
    const real_t vort = std::exp(-rad1 * oo_sigma2);

    return vort;
};

/**
 * @brief returns the value of a scalar compactly supported tube aligned with the normal direction
 * 
 */
lambda_t<real_t, const real_t[], const real_t[], const real_t, const real_t, const lda_t>
    scalar_compact_tube = [](const real_t pos[3], const real_t center[3], const real_t sigma, const real_t beta, const lda_t normal) -> real_t {
    const lda_t idx = (normal + 1) % 3;
    const lda_t idy = (normal + 2) % 3;

    const real_t gamma     = beta * sigma;
    const real_t oo_sigma2 = 1.0 / (sigma * sigma);
    const real_t oo_gamma2 = 1.0 / (gamma * gamma);
    // compute the gaussian
    const real_t rad_sq = pow(pos[idx] - center[idx], 2) + pow(pos[idy] - center[idy], 2);
    const real_t rho1   = rad_sq * oo_sigma2;
    const real_t rho2   = rad_sq * oo_gamma2;
    const real_t vort   = (rho2 < 1.0) ? (exp(-rho1 / (1.0 - rho2))) : 0.0;
    return vort;
};

/**
 * @brief returns the value of a scalar ring perpendicular to the normal direction
 * 
 */
lambda_t<real_t, const real_t[], const real_t[], const real_t, const real_t, const lda_t>
    scalar_ring = [](const real_t pos[3], const real_t center[3],
                     const real_t radius, const real_t sigma, const lda_t normal) -> real_t {
    const lda_t idx = (normal + 1) % 3;
    const lda_t idy = (normal + 2) % 3;
    const lda_t idz = normal;

    const real_t oo_sigma2 = 1.0 / (sigma * sigma);

    // compute the gaussian
    const real_t alpha   = 2.0 * M_PI + atan2(pos[idy] - center[idy], pos[idx] - center[idx]);
    const real_t radr    = sqrt(pow(pos[idx] - center[idx], 2) + pow(pos[idy] - center[idy], 2));
    const real_t rad1r   = radr - radius;
    const real_t rad2r   = radr + radius;
    const real_t rad1_sq = pow(rad1r, 2) + pow(pos[idz] - center[idz], 2);
    const real_t rad2_sq = pow(rad2r, 2) + pow(pos[idz] - center[idz], 2);
    const real_t vort    = (std::exp(-rad1_sq * oo_sigma2) + std::exp(-rad2_sq * oo_sigma2));

    return vort;
};

/**
 * @brief returns the value of a scalar compact ring perpendicular to the normal direction
 * 
 * we rely on the compact gaussian formula: exp( - (r/sigma)^2 / (1 - (r/gamma)^2) )
 * with:
 * - sigma = standard variation of the Guassian
 * - gamma = compactness of the gaussian
 * 
 * we define beta = gamma/sigma
 * 
 */
lambda_t<real_t, const real_t[], const real_t[], const lda_t, const real_t, const real_t, const real_t>
    scalar_compact_ring = [](const real_t pos[3], const real_t center[3],
                             const lda_t normal, const real_t radius, const real_t sigma, const real_t beta) -> real_t {
    //-------------------------------------------------------------------------
    m_assert(beta > 1, "the beta parameter = %f must be >1", beta);
    m_assert((beta * sigma) < radius, "to avoid cross combinatin, we must have that alpha < radius");

    const lda_t idx = (normal + 1) % 3;
    const lda_t idy = (normal + 2) % 3;
    const lda_t idz = normal;

    const real_t gamma     = beta * sigma;
    const real_t oo_sigma2 = 1.0 / (sigma * sigma);
    const real_t oo_gamma2 = 1.0 / (gamma * gamma);

    const real_t x_center = center[idx];
    const real_t y_center = center[idy];
    const real_t x        = pos[idx] - x_center;
    const real_t y        = pos[idy] - y_center;

    const real_t alpha = atan2(y, x);

    // add the modes
    real_t z_center = center[idz];
    // const short_t n_mode   = freq_rad.size();
    // for (short_t id = 0; id < n_mode; ++id) {
    //     z_center += amp_rad[id] * radius * 1.0 / n_mode * sin(freq_rad[id] * alpha);
    // }
    const real_t z = pos[idz] - z_center;

    // compute the gaussian

    const real_t radr   = sqrt(pow(x, 2) + pow(y, 2));
    const real_t rad    = radr - radius;
    const real_t rad_sq = pow(rad, 2) + pow(z, 2);
    const real_t rho1   = rad_sq * oo_sigma2;
    const real_t rho2   = rad_sq * oo_gamma2;
    const real_t vort   = (rho2 < 1.0) ? exp(-rho1 / (1.0 - rho2)) : 0.0;

    return vort;
    //-------------------------------------------------------------------------
};

// void SetValue::operator()(const ForestGrid*  grid, Field*  field, const lda_t ida) {
//     m_begin;
//     //-------------------------------------------------------------------------
//     // get the span of ida
//     ida_start_ = ida;
//     ida_end_   = ida + 1;
//     // go for it
//     DoOpMesh(this, &SetValue::FillGridBlock, grid, field);
//     // update the ghost status
//     m_verb("setting the ghosts of %s to false", field->name().c_str());
//     // we cannot set the ghost status as only one direction has been done...
//     field->ghost_status(false);
//     //-------------------------------------------------------------------------
//     m_end;
// }

// //=====================================================================================================
// SetAbs::SetAbs(const real_t alpha[3], const real_t center[3]) : SetAbs(alpha, center, nullptr) {}
// SetAbs::SetAbs(const real_t alpha[3], const real_t center[3], const Wavelet* interp) : SetValue(interp) {
//     m_begin;
//     //-------------------------------------------------------------------------
//     for (lda_t id = 0; id < 3; id++) {
//         center_[id] = center[id];
//         alpha_[id]  = alpha[id];
//     }
//     //-------------------------------------------------------------------------
//     m_end;
// }

// void SetAbs::FillGridBlock(const qid_t* qid, GridBlock* block, Field* fid) {
//     //-------------------------------------------------------------------------
//     real_t        pos[3];
//     const real_t* xyz   = block->xyz();
//     const real_t* hgrid = block->hgrid();

//     data_ptr block_data = block->data(fid);
//     for (lda_t ida = ida_start_; ida < ida_end_; ida++) {
//         real_t* data = block_data.Write(ida, block);
//         for (lid_t i2 = start_; i2 < end_; i2++) {
//             for (lid_t i1 = start_; i1 < end_; i1++) {
//                 for (lid_t i0 = start_; i0 < end_; i0++) {
//                     // get the position
//                     m_pos(pos, i0, i1, i2, hgrid, xyz);
//                     // compute the gaussian
//                     const real_t x0 = pos[0] - center_[0];
//                     const real_t x1 = pos[1] - center_[1];
//                     const real_t x2 = pos[2] - center_[2];

//                     data[m_idx(i0, i1, i2)] = alpha_[0] * std::fabs(x0) + alpha_[1] * std::fabs(x1) + alpha_[2] * std::fabs(x2);
//                 }
//             }
//         }
//     }
//     //-------------------------------------------------------------------------
// }

// //=====================================================================================================
// SetSinus::SetSinus(const real_t length[3], const real_t freq[3], const real_t alpha[3]) : SetSinus(length, freq, alpha, nullptr) {}
// SetSinus::SetSinus(const real_t length[3], const real_t freq[3], const real_t alpha[3], const Wavelet*  interp) : SetValue(interp) {
//     m_begin;
//     //-------------------------------------------------------------------------
//     for (lda_t id = 0; id < 3; id++) {
//         length_[id] = length[id];
//         freq_[id]   = freq[id];
//         alpha_[id]  = alpha[id];
//     }
//     //-------------------------------------------------------------------------
//     m_end;
// }

// void SetSinus::FillGridBlock(const qid_t*  qid, GridBlock*  block, Field*  fid) {
//     //-------------------------------------------------------------------------
//     real_t        pos[3];
//     const real_t* xyz     = block->xyz();
//     const real_t* hgrid   = block->hgrid();
//     const real_t  fact[3] = {2.0 * M_PI * freq_[0] / length_[0], 2.0 * M_PI * freq_[1] / length_[1], 2.0 * M_PI * freq_[2] / length_[2]};

//     data_ptr block_data = block->data(fid);
//     for (lda_t ida = ida_start_; ida < ida_end_; ida++) {
//         real_t* data = block_data.Write(ida, block);
//         for (lid_t i2 = start_; i2 < end_; i2++) {
//             for (lid_t i1 = start_; i1 < end_; i1++) {
//                 for (lid_t i0 = start_; i0 < end_; i0++) {
//                     // get the position
//                     m_pos(pos, i0, i1, i2, hgrid, xyz);
//                     data[m_idx(i0, i1, i2)] = alpha_[0] * sin(pos[0] * fact[0]) +
//                                               alpha_[1] * sin(pos[1] * fact[1]) +
//                                               alpha_[2] * sin(pos[2] * fact[2]);
//                 }
//             }
//         }
//     }
//     //-------------------------------------------------------------------------
// }

// //=====================================================================================================
// SetCosinus::SetCosinus(const real_t length[3], const real_t freq[3], const real_t alpha[3]) : SetCosinus(length, freq, alpha, nullptr) {}
// SetCosinus::SetCosinus(const real_t length[3], const real_t freq[3], const real_t alpha[3], const Wavelet*  interp) : SetValue(interp) {
//     m_begin;
//     //-------------------------------------------------------------------------
//     for (int id = 0; id < 3; id++) {
//         length_[id] = length[id];
//         freq_[id]   = freq[id];
//         alpha_[id]  = alpha[id];
//     }
//     //-------------------------------------------------------------------------
//     m_end;
// }

// void SetCosinus::FillGridBlock(const qid_t*  qid, GridBlock*  block, Field*  fid) {
//     //-------------------------------------------------------------------------
//     real_t        pos[3];
//     const real_t* xyz     = block->xyz();
//     const real_t* hgrid   = block->hgrid();
//     const real_t  vol     = 1.0 / (hgrid[0] * hgrid[1] * hgrid[2]);
//     const real_t  fact[3] = {2.0 * M_PI * freq_[0] / length_[0], 2.0 * M_PI * freq_[1] / length_[1], 2.0 * M_PI * freq_[2] / length_[2]};

//     data_ptr block_data = block->data(fid);
//     for (lda_t ida = ida_start_; ida < ida_end_; ida++) {
//         real_t* data = block_data.Write(ida, block);
//         for (lid_t i2 = start_; i2 < end_; i2++) {
//             for (lid_t i1 = start_; i1 < end_; i1++) {
//                 for (lid_t i0 = start_; i0 < end_; i0++) {
//                     // get the position
//                     m_pos(pos, i0, i1, i2, hgrid, xyz);
//                     data[m_idx(i0, i1, i2)] = alpha_[0] * cos(pos[0] * fact[0]) +
//                                               alpha_[1] * cos(pos[1] * fact[1]) +
//                                               alpha_[2] * cos(pos[2] * fact[2]);
//                 }
//             }
//         }
//     }
//     //-------------------------------------------------------------------------
// }

// //=====================================================================================================
// SetPolynom::SetPolynom(const lid_t degree[3], const real_t direction[3], const real_t shift[3]) : SetPolynom(degree, direction, shift, nullptr) {}
// SetPolynom::SetPolynom(const lid_t degree[3], const real_t direction[3], const real_t shift[3], const Wavelet*  interp) : SetValue(interp) {
//     m_begin;
//     //-------------------------------------------------------------------------
//     for (lda_t id = 0; id < 3; id++) {
//         deg_[id]   = degree[id];
//         dir_[id]   = direction[id];
//         shift_[id] = shift[id];
//     }
//     //-------------------------------------------------------------------------
//     m_end;
// }
// void SetPolynom::FillGridBlock(const qid_t*  qid, GridBlock*  block, Field*  fid) {
//     //-------------------------------------------------------------------------
//     real_t        pos[3];
//     const real_t* xyz   = block->xyz();
//     const real_t* hgrid = block->hgrid();

//     data_ptr block_data = block->data(fid);
//     for (lda_t ida = ida_start_; ida < ida_end_; ida++) {
//         real_t* data = block_data.Write(ida, block);
//         for (lid_t i2 = start_; i2 < end_; i2++) {
//             for (lid_t i1 = start_; i1 < end_; i1++) {
//                 for (lid_t i0 = start_; i0 < end_; i0++) {
//                     // get the position
//                     m_pos(pos, i0, i1, i2, hgrid, xyz);

//                     data[m_idx(i0, i1, i2)] = dir_[0] * pow(pos[0] - shift_[0], deg_[0]) +
//                                               dir_[1] * pow(pos[1] - shift_[1], deg_[1]) +
//                                               dir_[2] * pow(pos[2] - shift_[2], deg_[2]);
//                 }
//             }
//         }
//     }
//     //-------------------------------------------------------------------------
// }

// //=====================================================================================================
// SetExponential::SetExponential(const real_t center[3], const real_t sigma[3], const real_t alpha) : SetExponential(center, sigma, alpha, nullptr) {}
// SetExponential::SetExponential(const real_t center[3], const real_t sigma[3], const real_t alpha, const Wavelet*  interp) : SetValue(interp) {
//     m_begin;
//     //-------------------------------------------------------------------------
//     for (int id = 0; id < 3; id++) {
//         center_[id] = center[id];
//         sigma_[id]  = sigma[id];
//     }
//     alpha_ = alpha;
//     //-------------------------------------------------------------------------
//     m_end;
// }

// void SetExponential::FillGridBlock(const qid_t*  qid, GridBlock*  block, Field*  fid) {
//     //-------------------------------------------------------------------------
//     real_t        pos[3];
//     const real_t* xyz   = block->xyz();
//     const real_t* hgrid = block->hgrid();

//     real_t sigma     = sqrt(sigma_[0] * sigma_[0] + sigma_[1] * sigma_[1] + sigma_[2] * sigma_[2]);
//     real_t oo_sigma2 = (sigma > 0.0) ? 1.0 / (sigma * sigma) : 0.0;
//     real_t fact      = alpha_ * sqrt(1.0 / M_PI * oo_sigma2);

//     data_ptr block_data = block->data(fid);
//     for (lda_t ida = ida_start_; ida < ida_end_; ida++) {
//         real_t* data = block_data.Write(ida, block);
//         for (lid_t i2 = start_; i2 < end_; i2++) {
//             for (lid_t i1 = start_; i1 < end_; i1++) {
//                 for (lid_t i0 = start_; i0 < end_; i0++) {
//                     // get the position
//                     m_pos(pos, i0, i1, i2, hgrid, xyz);
//                     // compute the gaussian
//                     const real_t rhox = (sigma_[0] > 0) ? ((pos[0] - center_[0]) / sigma_[0]) : 0.0;
//                     const real_t rhoy = (sigma_[1] > 0) ? ((pos[1] - center_[1]) / sigma_[1]) : 0.0;
//                     const real_t rhoz = (sigma_[2] > 0) ? ((pos[2] - center_[2]) / sigma_[2]) : 0.0;
//                     const real_t rho  = rhox * rhox + rhoy * rhoy + rhoz * rhoz;

//                     data[m_idx(i0, i1, i2)] = fact * std::exp(-rho);
//                 }
//             }
//         }
//     }
//     //-------------------------------------------------------------------------
// }

// //=====================================================================================================
// SetErf::SetErf(const real_t center[3], const real_t sigma[3], const real_t alpha) : SetErf(center, sigma, alpha, nullptr) {}
// SetErf::SetErf(const real_t center[3], const real_t sigma[3], const real_t alpha, const Wavelet*  interp) : SetValue(interp) {
//     m_begin;
//     //-------------------------------------------------------------------------
//     for (lda_t id = 0; id < 3; id++) {
//         center_[id] = center[id];
//         sigma_[id]  = sigma[id];
//     }
//     alpha_ = alpha;
//     //-------------------------------------------------------------------------
//     m_end;
// }

// void SetErf::FillGridBlock(const qid_t*  qid, GridBlock*  block, Field*  fid) {
//     //-------------------------------------------------------------------------
//     real_t        pos[3];
//     const real_t* xyz   = block->xyz();
//     const real_t* hgrid = block->hgrid();

//     real_t sigma     = sqrt(sigma_[0] * sigma_[0] + sigma_[1] * sigma_[1] + sigma_[2] * sigma_[2]);
//     real_t oo_sigma2 = (sigma > 0.0) ? 1.0 / (sigma * sigma) : 0.0;
//     real_t oo_sqrt2  = 1.0 / M_SQRT2;
//     real_t fact      = alpha_ / (4.0 * M_PI * sigma);  // see Wincky encyclopedia

//     data_ptr block_data = block->data(fid);
//     for (lda_t ida = ida_start_; ida < ida_end_; ida++) {
//         real_t* data = block_data.Write(ida, block);
//         for (lid_t i2 = start_; i2 < end_; i2++) {
//             for (lid_t i1 = start_; i1 < end_; i1++) {
//                 for (lid_t i0 = start_; i0 < end_; i0++) {
//                     // get the position
//                     m_pos(pos, i0, i1, i2, hgrid, xyz);
//                     // compute the gaussian
//                     const real_t x0         = (sigma_[0] > 0) ? pos[0] - center_[0] : 0.0;
//                     const real_t x1         = (sigma_[1] > 0) ? pos[1] - center_[1] : 0.0;
//                     const real_t x2         = (sigma_[2] > 0) ? pos[2] - center_[2] : 0.0;
//                     const real_t rho2       = (x0 * x0 + x1 * x1 + x2 * x2) * oo_sigma2;
//                     const real_t rho        = sqrt(rho2);
//                     data[m_idx(i0, i1, i2)] = fact / rho * std::erf(rho * oo_sqrt2);
//                 }
//             }
//         }
//     }
//     //-------------------------------------------------------------------------
// }

// //=====================================================================================================
// SetVortexRing::SetVortexRing(const lda_t normal, const real_t center[3], const real_t sigma, const real_t radius) : SetVortexRing(normal, center, sigma, radius, nullptr) {}
// SetVortexRing::SetVortexRing(const lda_t normal, const real_t center[3], const real_t sigma, const real_t radius, const Wavelet*  interp) : SetValue(interp) {
//     m_begin;
//     //-------------------------------------------------------------------------
//     normal_ = normal;
//     sigma_  = sigma;
//     radius_ = radius;
//     for (lda_t id = 0; id < 3; ++id) {
//         center_[id] = center[id];
//     }
//     //-------------------------------------------------------------------------
//     m_end;
// }

// void SetVortexRing::FillGridBlock(const qid_t*  qid, GridBlock*  block, Field*  fid) {
//     //-------------------------------------------------------------------------
//     real_t        pos[3];
//     const real_t* xyz   = block->xyz();
//     const real_t* hgrid = block->hgrid();

//     const real_t oo_sigma2   = 1.0 / (sigma_ * sigma_);
//     const real_t oo_pisigma2 = 1.0;  /// sqrt(M_PI * sigma_ * sigma_); //todo change that because sqrt(M_PI * sigma_ * sigma_) is the initial amplitude

//     // compute the normal direction as being the z one and the two other as x and y
//     const lda_t idx = (normal_ + 1) % 3;
//     const lda_t idy = (normal_ + 2) % 3;
//     const lda_t idz = normal_;
//     // get the pointers correct
//     data_ptr block_data = block->data(fid);
//     real_t*  wx         = block_data.Write(idx, block);
//     real_t*  wy         = block_data.Write(idy, block);
//     real_t*  wz         = block_data.Write(idz, block);

//     for (lid_t i2 = start_; i2 < end_; i2++) {
//         for (lid_t i1 = start_; i1 < end_; i1++) {
//             for (lid_t i0 = start_; i0 < end_; i0++) {
//                 // get the position
//                 m_pos(pos, i0, i1, i2, hgrid, xyz);
//                 // pos[0] = fmod(pos[0]+1, 1.0);
//                 // pos[1] = fmod(pos[1]+1, 1.0);
//                 // pos[2] = fmod(pos[2]+2, 1.0);

//                 // wrt to the center
//                 const real_t alpha   = 2.0 * M_PI + atan2(pos[idy] - center_[idy], pos[idx] - center_[idx]);
//                 const real_t radr    = sqrt(pow(pos[idx] - center_[idx], 2) + pow(pos[idy] - center_[idy], 2));
//                 const real_t rad1r   = radr - radius_;
//                 const real_t rad2r   = radr + radius_;
//                 const real_t rad1_sq = pow(rad1r, 2) + pow(pos[idz] - center_[idz], 2);
//                 const real_t rad2_sq = pow(rad2r, 2) + pow(pos[idz] - center_[idz], 2);
//                 const real_t vort    = oo_pisigma2 * (exp(-rad1_sq * oo_sigma2) - exp(-rad2_sq * oo_sigma2));

//                 wx[m_idx(i0, i1, i2)] = -vort * sin(alpha);  //  -vort * cos(alpha);
//                 wy[m_idx(i0, i1, i2)] = vort * cos(alpha);   // +vort * sin(alpha);
//                 wz[m_idx(i0, i1, i2)] = 0.0;

//                 // // wrt to the center
//                 // const real_t alpha = atan2(pos[idy] - center_[idy], pos[idx] - center_[idx]);
//                 // const real_t x     = pos[idx] - (center_[idx] + radius_ * cos(alpha));
//                 // const real_t y     = pos[idy] - (center_[idy] + radius_ * sin(alpha));
//                 // const real_t z     = pos[idz] - (center_[idz]);
//                 // // get the local coords

//                 // // const real_t r_plane = sqrt(x * x + y * y) - radius_;
//                 // const real_t r2   = (x * x + y * y + z * z);
//                 // const real_t rho2 = r2 * oo_sigma2;
//                 // const real_t vort = oo_pisigma2 * exp(-rho2);
//                 // // const real_t vort = oo_pisigma2 * exp(1.0 - 1.0 / (1.0 - r2 * oo_sigma2));
//                 // wx[m_idx(i0, i1, i2)] = -vort * sin(alpha);  //  -vort * cos(alpha);
//                 // wy[m_idx(i0, i1, i2)] = vort * cos(alpha);   // +vort * sin(alpha);
//                 // wz[m_idx(i0, i1, i2)] = 0.0;
//             }
//         }
//     }
//     //-------------------------------------------------------------------------
// }

// //=====================================================================================================
// SetCompactVortexRing::SetCompactVortexRing(const lda_t normal, const real_t center[3], const real_t sigma, const real_t radius, const real_t cutoff) : SetCompactVortexRing(normal, center, sigma, radius, cutoff, nullptr) {}
// SetCompactVortexRing::SetCompactVortexRing(const lda_t normal, const real_t center[3], const real_t sigma, const real_t radius, const real_t cutoff, const Wavelet*  interp) : SetValue(interp) {
//     m_begin;
//     //-------------------------------------------------------------------------
//     normal_ = normal;
//     sigma_  = sigma;
//     radius_ = radius;
//     cutoff_ = cutoff;
//     for (lda_t id = 0; id < 3; ++id) {
//         center_[id] = center[id];
//     }
//     //-------------------------------------------------------------------------
//     m_end;
// }

// void SetCompactVortexRing::FillGridBlock(const qid_t*  qid, GridBlock*  block, Field*  fid) {
//     //-------------------------------------------------------------------------
//     real_t        pos[3];
//     const real_t* xyz   = block->xyz();
//     const real_t* hgrid = block->hgrid();

//     const real_t oo_sigma2   = 1.0 / (sigma_ * sigma_);
//     const real_t oo_R2       = 1.0 / (cutoff_ * cutoff_);
//     const real_t oo_pisigma2 = 1.0 / (M_PI * sigma_ * sigma_);

//     // compute the normal direction as being the z one and the two other as x and y
//     const lda_t idx = (normal_ + 1) % 3;
//     const lda_t idy = (normal_ + 2) % 3;
//     const lda_t idz = normal_;
//     // get the pointers correct
//     data_ptr block_data = block->data(fid);
//     real_t*  wx         = block_data.Write(idx, block);
//     real_t*  wy         = block_data.Write(idy, block);
//     real_t*  wz         = block_data.Write(idz, block);

//     for (lid_t i2 = start_; i2 < end_; i2++) {
//         for (lid_t i1 = start_; i1 < end_; i1++) {
//             for (lid_t i0 = start_; i0 < end_; i0++) {
//                 // get the position
//                 m_pos(pos, i0, i1, i2, hgrid, xyz);
//                 // wrt to the center
//                 const real_t alpha = atan2(pos[idy] - center_[idy], pos[idx] - center_[idx]);
//                 const real_t x     = pos[idx] - (center_[idx] + radius_ * cos(alpha));
//                 const real_t y     = pos[idy] - (center_[idy] + radius_ * sin(alpha));
//                 const real_t z     = pos[idz] - (center_[idz]);
//                 // get the local coords

//                 // const real_t r_plane = sqrt(x * x + y * y) - radius_;
//                 const real_t r2    = (x * x + y * y + z * z);
//                 const real_t rho2  = r2 * oo_sigma2;
//                 const real_t rhoR2 = r2 * oo_R2;

//                 if (rhoR2 < 1.0) {
//                     const real_t vort = oo_pisigma2 * exp(-rho2 / (1.0 - rhoR2));
//                     // const real_t vort = oo_pisigma2 * exp(1.0 - 1.0 / (1.0 - r2 * oo_sigma2));
//                     wx[m_idx(i0, i1, i2)] = -vort * sin(alpha);  //  -vort * cos(alpha);
//                     wy[m_idx(i0, i1, i2)] = vort * cos(alpha);   // +vort * sin(alpha);
//                     wz[m_idx(i0, i1, i2)] = 0.0;

//                 } else {
//                     wx[m_idx(i0, i1, i2)] = 0.0;
//                     wy[m_idx(i0, i1, i2)] = 0.0;
//                     wz[m_idx(i0, i1, i2)] = 0.0;
//                 }
//             }
//         }
//     }
//     //-------------------------------------------------------------------------
// }

// //=====================================================================================================
// SetABSVelocity::SetABSVelocity(const real_t a, const real_t b, const real_t c) : SetABSVelocity(a, b, c, nullptr) {}
// SetABSVelocity::SetABSVelocity(const real_t a, const real_t b, const real_t c, const Wavelet*  interp) : SetValue(interp) {
//     m_begin;
//     //-------------------------------------------------------------------------
//     a_ = a;
//     b_ = b;
//     c_ = c;
//     //-------------------------------------------------------------------------
//     m_end;
// }

// void SetABSVelocity::FillGridBlock(const qid_t*  qid, GridBlock*  block, Field*  fid) {
//     //-------------------------------------------------------------------------
//     real_t* vx = block->data(fid, 0).Write();
//     real_t* vy = block->data(fid, 1).Write();
//     real_t* vz = block->data(fid, 2).Write();

//     auto op = [=, &vx, &vy, &vz](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
//         real_t pos[3];
//         m_pos(pos, i0, i1, i2, block->hgrid(), block->xyz());

//         // get the params
//         constexpr real_t two_pi = 2.0 * M_PI;

//         // get the domain-relative position
//         const real_t x = pos[0];
//         const real_t y = pos[1];
//         const real_t z = pos[2];

//         vx[m_idx(i0, i1, i2)] = a_ * sin(two_pi * z) + c_ * cos(two_pi * y);
//         vy[m_idx(i0, i1, i2)] = b_ * sin(two_pi * x) + a_ * cos(two_pi * z);
//         vz[m_idx(i0, i1, i2)] = c_ * sin(two_pi * y) + b_ * cos(two_pi * x);
//     };

//     for_loop(&op, start_, end_);
//     //-------------------------------------------------------------------------
// }

// //=====================================================================================================
// SetScalarRing::SetScalarRing(const lda_t normal, const real_t center[3], const real_t sigma, const real_t radius, const real_t vel[3]) : SetScalarRing(normal, center, sigma, radius, vel, nullptr) {}
// SetScalarRing::SetScalarRing(const lda_t normal, const real_t center[3], const real_t sigma, const real_t radius, const real_t vel[3], const Wavelet*  interp) : SetValue(interp) {
//     m_begin;
//     //-------------------------------------------------------------------------
//     normal_ = normal;
//     sigma_  = sigma;
//     radius_ = radius;
//     for (lda_t id = 0; id < 3; ++id) {
//         center_[id] = center[id];
//         vel_[id] = vel[id];
//     }
//     //-------------------------------------------------------------------------
//     m_end;
// }

// void SetScalarRing::FillGridBlock(const qid_t*  qid, GridBlock*  block, Field*  fid) {
//     m_assert(fid->lda() == 1, "this function is for scalar fields only");
//     //-------------------------------------------------------------------------
//     real_t pos[3];

//     const real_t oo_sigma2   = 1.0 / (sigma_ * sigma_);
//     const real_t oo_pisigma2 = 1.0;  /// sqrt(M_PI * sigma_ * sigma_); //todo change that because sqrt(M_PI * sigma_ * sigma_) is the initial amplitude

//     // compute the normal direction as being the z one and the two other as x and y
//     const lda_t idx = (normal_ + 1) % 3;
//     const lda_t idy = (normal_ + 2) % 3;
//     const lda_t idz = normal_;

//     const real_t center[3] = {center_[0] + time_ * vel_[0],
//                               center_[1] + time_ * vel_[1],
//                               center_[2] + time_ * vel_[2]};

//     // get the pointers correct
//     real_t* w = block->data(fid).Write();

//     auto op = [=, &w](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
//         // get the position
//         real_t pos[3];
//         m_pos(pos, i0, i1, i2, block->hgrid(), block->xyz());

//         // wrt to the center
//         const real_t alpha   = 2.0 * M_PI + atan2(pos[idy] - center[idy], pos[idx] - center[idx]);
//         const real_t radr    = sqrt(pow(pos[idx] - center[idx], 2) + pow(pos[idy] - center[idy], 2));
//         const real_t rad1r   = radr - radius_;
//         const real_t rad2r   = radr + radius_;
//         const real_t rad1_sq = pow(rad1r, 2) + pow(pos[idz] - center[idz], 2);
//         const real_t rad2_sq = pow(rad2r, 2) + pow(pos[idz] - center[idz], 2);
//         const real_t vort    = oo_pisigma2 * (exp(-rad1_sq * oo_sigma2) + exp(-rad2_sq * oo_sigma2));

//         w[m_idx(i0, i1, i2)] = vort;
//     };

//     for_loop(&op, start_, end_);
//     //-------------------------------------------------------------------------
// }

// //=====================================================================================================
// SetScalarTube::SetScalarTube(const lda_t dir, const real_t center[3], const real_t sigma, const real_t b) : SetScalarTube(dir, center, sigma, b, nullptr) {}
// SetScalarTube::SetScalarTube(const lda_t dir, const real_t center[3], const real_t sigma, const real_t b, const Wavelet*  interp) : SetValue(interp) {
//     m_begin;
//     //-------------------------------------------------------------------------
//     dir_   = dir;
//     sigma_ = sigma;
//     b_     = b;
//     for (lda_t id = 0; id < 3; ++id) {
//         center_[id] = center[id];
//     }
//     //-------------------------------------------------------------------------
//     m_end;
// }

// void SetScalarTube::FillGridBlock(const qid_t*  qid, GridBlock*  block, Field*  fid) {
//     m_assert(fid->lda() == 1, "this function is for scalar fields only");
//     //-------------------------------------------------------------------------
//     real_t        pos[3];
//     const real_t* xyz   = block->xyz();
//     const real_t* hgrid = block->hgrid();

//     const real_t oo_sigma2   = 1.0 / (sigma_ * sigma_);
//     const real_t oo_pisigma2 = 1.0;  /// sqrt(M_PI * sigma_ * sigma_); //todo change that because sqrt(M_PI * sigma_ * sigma_) is the initial amplitude

//     // compute the normal direction as being the z one and the two other as x and y
//     const lda_t idx = (dir_ + 1) % 3;
//     const lda_t idy = (dir_ + 2) % 3;
//     const lda_t idz = dir_;

//     // get the pointers correct
//     real_t* w0 = block->data(fid, 0).Write();
//     real_t* w1 = block->data(fid, 1).Write();

//     auto op = [=, &w0, &w1](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
//         // get the position
//         real_t pos[3];
//         m_pos(pos, i0, i1, i2, block->hgrid(), block->xyz());

//         // wrt to the center
//         // const real_t alpha   = 2.0 * M_PI + atan2(pos[idy] - center_[idy], pos[idx] - center_[idx]);
//         // const real_t radr    = sqrt(pow(pos[idx] - center_[idx], 2) + pow(pos[idy] - center_[idy], 2));
//         // const real_t rad1 = pow(pos[idx] - (center_[idx] + b_), 2) + pow(pos[idy] - center_[idy], 2);
//         // const real_t rad2 = pow(pos[idx] - (center_[idx] - b_), 2) + pow(pos[idy] - center_[idy], 2);
//         const real_t rad1 = pow(pos[idx] - (center_[idx]), 2) + pow(pos[idy] - center_[idy], 2);

//         // const real_t rad1_sq = pow(rad1r, 2) + pow(pos[idz] - center_[idz], 2);
//         // const real_t rad2_sq = pow(rad2r, 2) + pow(pos[idz] - center_[idz], 2);
//         // const real_t vort     = oo_pisigma2 * (exp(-rad1 * oo_sigma2) - exp(-rad2 * oo_sigma2));
//         w0[m_idx(i0, i1, i2)] = m_max(0.0, oo_pisigma2 * exp(-rad1 * oo_sigma2));
//         // w1[m_idx(i0, i1, i2)] = m_max(0.0, oo_pisigma2 * exp(-rad2 * oo_sigma2));
//         // w2[m_idx(i0, i1, i2)] = m_max(0.0, oo_pisigma2 * exp(-rad2 * oo_sigma2));
//     };

//     for_loop(&op, start_, end_);
//     //-------------------------------------------------------------------------
// }