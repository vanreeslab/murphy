#include "advection.hpp"


/**
 * @brief Apply the WENO scheme of the third order (r=2, 2*r-1 = 3)
 * 
 * See Jiang1996, Schu1997, and Don2013 for more details
 * From Shu1997 p4 we know that
 * ```
 *      -1/2 3/2
 *            1/2 1/2
 * ```
 * 
 * @tparam  
 * @param qid 
 * @param block 
 * @param is_outer 
 * @param fid_src 
 * @param fid_trg 
 */
template <>
void Advection<M_WENO_Z, 3>::DoMagic(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, const bool is_outer, m_ptr<const Field> fid_src, m_ptr<Field> fid_trg) const {
    m_assert(u_->ghost_status(), "the ghost values of the velocity must be known!");
    //-------------------------------------------------------------------------
    const real_t h[3] = {block->hgrid(0), block->hgrid(1), block->hgrid(2)};
    // get the lambda to compute a flux: ida = direction of the flux, fvel = velocity on the face
    // the flux computed is the one in i+1/2!
    auto flux_weno_3 = [=](const lda_t ida, const real_t* src, const real_t* vel) -> real_t {
        // get the face velocity
        const real_t fvel = 0.5 * (vel[m_idx_delta(+1, ida)] + vel[0]);
        // get the values
        const real_t fm1 = src[m_idx_delta(-1, ida)] * vel[m_idx_delta(-1, ida)];
        const real_t fm0 = src[0] * vel[0];
        const real_t fp1 = src[m_idx_delta(+1, ida)] * vel[m_idx_delta(+1, ida)];
        const real_t fp2 = src[m_idx_delta(+2, ida)] * vel[m_idx_delta(+2, ida)];

        // stencil 0 -> r = 1
        const real_t gamma_0 = 1.0 / 3.0;
        const real_t beta_0  = pow(fm0 - fm1, 2);
        const real_t s0      = (-0.5 * fm1 + 1.5 * fm0);
        // stencil 1 -> r = 0
        const real_t gamma_1 = 2.0 / 3.0;
        const real_t beta_1  = pow(fp1 - fm0, 2);
        const real_t s1      = (0.5 * fm0 + 0.5 * fp1);
        // stencil 2 -> r = -1
        const real_t gamma_2 = 1.0 / 3.0;
        const real_t beta_2  = pow(fp2 - fp1, 2);
        const real_t s2      = (1.5 * fp1 - 0.5 * fp2);

        // case 1: face velocity is > 0, we use stencil 0 and stencil 1
        real_t flux = 0.0;
        {
            const real_t tau     = fabs(beta_0 - beta_1);
            const real_t alpha_0 = gamma_0 * (1.0 + pow(tau / (beta_0 + pow(h[ida], 4)), 2));
            const real_t alpha_1 = gamma_1 * (1.0 + pow(tau / (beta_1 + pow(h[ida], 4)), 2));
            const real_t denom   = 1.0 / (alpha_0 + alpha_1);
            const real_t w0      = alpha_0 * denom;
            const real_t w1      = alpha_1 * denom;
            flux                 = flux + (m_sign(fvel) > 0.0) * (w0 * s0 + w1 * s1);
        }
        // case 2: face velocity is < 0, we use stencil 1 and 2
        {
            const real_t tau_3   = fabs(beta_1 - beta_2);
            const real_t alpha_1 = gamma_1 * (1.0 + pow(tau_3 / (beta_1 + pow(h[ida], 4)), 2));
            const real_t alpha_2 = gamma_2 * (1.0 + pow(tau_3 / (beta_2 + pow(h[ida], 4)), 2));
            const real_t denom   = 1.0 / (alpha_1 + alpha_2);
            const real_t w1      = alpha_1 * denom;
            const real_t w2      = alpha_2 * denom;
            flux                 = flux + (m_sign(fvel) < 0.0) * (w1 * s1 + w2 * s2);
        }
        return flux;
    };
    //-------------------------------------------------------------------------
    // get the data
    const real_t* data_u   = block->data(u_, 0).Read();
    const real_t* data_v   = block->data(u_, 1).Read();
    const real_t* data_w   = block->data(u_, 2).Read();
    const real_t* data_src = block->data(fid_src, ida_).Read();
    real_t*       data_trg = block->data(fid_trg, ida_).Write();

    // usefull factors
    const real_t alpha    = (accumulate_) ? 1.0 : 0.0;
    const real_t oneoh[3] = {1.0 / block->hgrid(0), 1.0 / block->hgrid(1), 1.0 / block->hgrid(2)};

    auto op = [=, &data_trg](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
        const bidx_t  idx[3] = {i0, i1, i2};
        real_t*       trg    = data_trg + m_idx(i0, i1, i2);
        const real_t* src    = data_src + m_idx(i0, i1, i2);
        const real_t* vel[3] = {data_u + m_idx(i0, i1, i2),
                                data_v + m_idx(i0, i1, i2),
                                data_w + m_idx(i0, i1, i2)};

        // reset the value if needed
        trg[0] *= alpha;

        // loop on the faces
        for (lda_t ida = 0; ida < 3; ++ida) {
            // as we are on the face +1/2, apply the flux on the:
            // -> left only if we are >=0
            // -> right only if we are < (M_N-1)
            const real_t flux_apply_left  = idx[ida] >= 0;
            const real_t flux_apply_right = idx[ida] < (M_N - 1);
            const real_t flux             = flux_weno_3(ida, src, vel[ida]);

            trg[m_idx_delta(0, ida)] -= oneoh[0] * flux * flux_apply_left;
            trg[m_idx_delta(+1, ida)] += oneoh[0] * flux * flux_apply_right;
        }
    };

    if (!is_outer) {
        for_loop<M_GS, M_N - M_GS>(&op);
    } else {
        // do the most on the X side to use vectorization
        for_loop<0, M_N, 0, M_N, 0, M_GS>(&op);          // Z-
        for_loop<0, M_N, 0, M_N, M_N - M_GS, M_N>(&op);  // Z+

        for_loop<0, M_N, 0, M_GS, M_GS, M_N - M_GS>(&op);          // Y-
        for_loop<0, M_N, M_N - M_GS, M_N, M_GS, M_N - M_GS>(&op);  // Y+

        for_loop<0, M_GS, M_GS, M_N - M_GS, M_GS, M_N - M_GS>(&op);          // X-
        for_loop<M_N - M_GS, M_N, M_GS, M_N - M_GS, M_GS, M_N - M_GS>(&op);  // X+
    }
    // -------------------------------------------------------------------------
}
