#include "advection.hpp"

#include <functional>

/**
 * @brief WENO-Z smoothing coefficient computation
 */
constexpr real_t alpha_z(const real_t gamma, const real_t beta, const real_t tau, const real_t h) {
    return gamma * (1.0 + pow(tau / (beta + pow(h, 4)), 2));
}

/**
 * @brief WENO-JS smoothing coefficient computation
 */
constexpr real_t alpha_js(const real_t gamma, const real_t beta) {
    constexpr real_t eps = 1e-8;
    return gamma * pow(1.0 / (eps + beta), 2);
}

/**
 * @brief defines a lambda that returns two fluxes: f- = flux[0] at (i-1/2, right biased) and f+ = flux[1] in (i+1/2, left biased)
 * 
 * @warning the notation of f+ and f- is ambiguous (compared to Shu1997), here we define f(v) = f-(v) + f+(v) where df-/dv < 0
 * and df+/dv > 0. In the case of the advection, we have f(v) = u v, so the separation is made based on the velocity u.
 * 
 */
using lambda_flux_weno_z_t = std::function<void(const lda_t ida, const real_t* const src, const real_t* const vel, real_t* const flux)>;

/**
 * @brief Given a flux function, do the magic for the WENO_Z schemes
 * 
 */
static void DoMagic_WENOZ(/* flux weno */ const lambda_flux_weno_z_t* flux_weno,
                          /* fields */ const real_t* data_vel[3], const real_t* data_src, real_t* data_trg,
                          /* factors */ const bool is_outer, const real_t alpha, const real_t oneoh[3]) {
    //-------------------------------------------------------------------------
    auto op = [=, &data_trg](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
        const bidx_t  idx[3] = {i0, i1, i2};
        real_t*       trg    = data_trg + m_idx(i0, i1, i2);
        const real_t* src    = data_src + m_idx(i0, i1, i2);

        // loop on the faces
        for (lda_t ida = 0; ida < 3; ++ida) {
            // as we are on the face +1/2, apply the flux on the:
            // -> left only if we are >0
            // -> right only if we are < (M_N-1)
            const bool flux_apply_left  = idx[ida] > 0;
            const bool flux_apply_right = idx[ida] < (M_N - 1);
            real_t     fluxes[2]        = {0.0, 0.0};
            (*flux_weno)(ida, src, data_vel[ida] + m_idx(i0, i1, i2), fluxes);

            // update
            trg[m_idx_delta(-1, ida)] -= oneoh[ida] * fluxes[0] * flux_apply_left;
            trg[m_idx_delta(+0, ida)] += oneoh[ida] * (fluxes[0] - fluxes[1]);
            trg[m_idx_delta(+1, ida)] += oneoh[ida] * fluxes[1] * flux_apply_right;
        }
    };

    auto reset = [=, &data_trg](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
        data_trg[m_idx(i0, i1, i2)] *= alpha;
    };

    if (!is_outer) {
        // reset the whole block, we don't need the ghost points to do that
        for_loop<0, M_N>(&reset);
        for_loop<M_GS, M_N - M_GS>(&op);
    } else {
        // do the most on the X side to use vectorization
        // need to start in -1 to put a flux in the point 0!
        for_loop<-1, M_N + 1, -1, M_N + 1, -1, M_GS>(&op);             // Z-
        for_loop<-1, M_N + 1, -1, M_N + 1, M_N - M_GS, M_N + 1>(&op);  // Z+

        for_loop<-1, M_N + 1, -1, M_GS, M_GS, M_N - M_GS>(&op);             // Y-
        for_loop<-1, M_N + 1, M_N - M_GS, M_N + 1, M_GS, M_N - M_GS>(&op);  // Y+

        for_loop<-1, M_GS, M_GS, M_N - M_GS, M_GS, M_N - M_GS>(&op);             // X-
        for_loop<M_N - M_GS, M_N + 1, M_GS, M_N - M_GS, M_GS, M_N - M_GS>(&op);  // X+
    }
    // -------------------------------------------------------------------------
}

/**
 * @brief Apply the WENO scheme of the third order (r=2, 2*r-1 = 3)
 * 
 * See Jiang1996, Schu1997, and Don2013 for more details
 * 
 * @tparam  
 * @param qid 
 * @param block 
 * @param is_outer 
 * @param fid_src 
 * @param fid_trg 
 */
template <>
void Advection<M_WENO_Z, 3>::DoMagic(const qid_t*  qid, GridBlock*  block, const bool is_outer, const Field*  fid_src, Field*  fid_trg) const {
    m_assert(u_->ghost_status(), "the ghost values of the velocity must be known!");
    //-------------------------------------------------------------------------
    const real_t h[3] = {block->hgrid(0), block->hgrid(1), block->hgrid(2)};

    // get the lambda to compute a flux: ida = direction of the flux, fvel = velocity on the face
    // the flux computed is the one in i+1/2!
    lambda_flux_weno_z_t flux_weno_3 = [=](const lda_t ida, const real_t* const src, const real_t* const vel, real_t* const flux) -> void {
        // get the values
        const real_t field[3] = {src[m_idx_delta(-1, ida)] * vel[m_idx_delta(-1, ida)],
                                 src[0] * vel[0],
                                 src[m_idx_delta(+1, ida)] * vel[m_idx_delta(+1, ida)]};
        // get the field shifted (it's easier)
        const real_t* const f = field + 1;

        // we compute the betas unique, they adapt to the stencils
        const real_t beta_0 = pow(f[-1] - f[0], 2);
        const real_t beta_1 = pow(f[0] - f[1], 2);
        const real_t tau_3  = fabs(beta_1 - beta_0);

        // case 0: inteface i-1/2, right biased stencils - flux is negative
        {
            const real_t gamma_0 = 2.0 / 3.0;
            const real_t s0      = (0.5 * f[-1] + 0.5 * f[0]);
            const real_t gamma_1 = 1.0 / 3.0;
            const real_t s1      = (1.5 * f[0] - 0.5 * f[1]);

            const real_t alpha_0 = alpha_z(gamma_0, beta_0, tau_3, h[ida]);
            const real_t alpha_1 = alpha_z(gamma_1, beta_1, tau_3, h[ida]);
            const real_t denom   = 1.0 / (alpha_0 + alpha_1);
            const real_t w0      = alpha_0 * denom;
            const real_t w1      = alpha_1 * denom;

            const real_t fvel = 0.5 * (vel[m_idx_delta(-1, ida)] + vel[0]);  // i - 1/2
            flux[0]           = (m_sign(fvel) < 0.0) * (w0 * s0 + w1 * s1);  // it's a negative flux
        }
        // case 1: inteface i+1/2, left biased stencils - flux is positive
        {
            const real_t gamma_0 = 1.0 / 3.0;
            const real_t s0      = (-0.5 * f[-1] + 1.5 * f[0]);
            const real_t gamma_1 = 2.0 / 3.0;
            const real_t s1      = (0.5 * f[0] + 0.5 * f[1]);

            const real_t alpha_0 = alpha_z(gamma_0, beta_0, tau_3, h[ida]);
            const real_t alpha_1 = alpha_z(gamma_1, beta_1, tau_3, h[ida]);
            const real_t denom   = 1.0 / (alpha_0 + alpha_1);
            const real_t w0      = alpha_0 * denom;
            const real_t w1      = alpha_1 * denom;

            const real_t fvel = 0.5 * (vel[m_idx_delta(+1, ida)] + vel[0]);  // i + 1/2
            flux[1]           = (m_sign(fvel) > 0.0) * (w0 * s0 + w1 * s1);  // it's a positive flux
        }
        // m_log("fluxes are %e %e",flux[0],flux[1]);
    };
    //-------------------------------------------------------------------------
    // apply the flux computations
    const real_t  alpha       = (accumulate_) ? 1.0 : 0.0;
    const real_t  oneoh[3]    = {1.0 / block->hgrid(0), 1.0 / block->hgrid(1), 1.0 / block->hgrid(2)};
    const real_t* data_vel[3] = {block->data(u_, 0).Read(), block->data(u_, 1).Read(), block->data(u_, 2).Read()};
    DoMagic_WENOZ(&flux_weno_3,
                  data_vel, block->data(fid_src, ida_).Read(), block->data(fid_trg, ida_).Write(),
                  is_outer, alpha, oneoh);
    // -------------------------------------------------------------------------
}

template <>
void Advection<M_WENO_Z, 5>::DoMagic(const qid_t*  qid, GridBlock*  block, const bool is_outer, const Field*  fid_src, Field*  fid_trg) const {
    m_assert(u_->ghost_status(), "the ghost values of the velocity must be known!");
    //-------------------------------------------------------------------------
    const real_t h[3] = {block->hgrid(0), block->hgrid(1), block->hgrid(2)};

    // get the lambda to compute a flux: ida = direction of the flux, fvel = velocity on the face
    lambda_flux_weno_z_t flux_weno_5 = [=](const lda_t ida, const real_t* src, const real_t* vel, real_t* const flux) -> void {
        // get the values
        const real_t  field[5] = {src[m_idx_delta(-2, ida)] * vel[m_idx_delta(-2, ida)],
                                 src[m_idx_delta(-1, ida)] * vel[m_idx_delta(-1, ida)],
                                 src[0] * vel[0],
                                 src[m_idx_delta(+1, ida)] * vel[m_idx_delta(+1, ida)],
                                 src[m_idx_delta(+2, ida)] * vel[m_idx_delta(+2, ida)]};
        const real_t* f        = field + 2;

        // we compute the betas unique, they adapt to the stencils
        const real_t beta_0 = 1.0 / 4.0 * pow(f[-2] - 4.0 * f[-1] + 3.0 * f[0], 2) + 13.0 / 12.0 * pow(f[-2] - 2.0 * f[-1] + f[0], 2);
        const real_t beta_1 = 1.0 / 4.0 * pow(-f[-1] + f[1], 2) + 13.0 / 12.0 * pow(f[-1] - 2.0 * f[0] + f[1], 2);
        const real_t beta_2 = 1.0 / 4.0 * pow(-3.0 * f[0] + 4.0 * f[1] - f[2], 2) + 13.0 / 12.0 * pow(f[0] - 2.0 * f[1] + f[2], 2);
        const real_t tau_5  = fabs(beta_2 - beta_0);

        // case 0: inteface i-1/2, right biased stencils - flux is negative
        {
            const real_t gamma_0 = 3.0 / 10.0;
            const real_t s0      = (-1.0 / 6.0 * f[-2] + 5.0 / 6.0 * f[-1] + 1.0 / 3.0 * f[0]);
            const real_t gamma_1 = 3.0 / 5.0;
            const real_t s1      = (1.0 / 3.0 * f[-1] + 5.0 / 6.0 * f[0] - 1.0 / 6.0 * f[1]);
            const real_t gamma_2 = 1.0 / 10.0;
            const real_t s2      = (11.0 / 6.0 * f[0] - 7.0 / 6.0 * f[1] + 1.0 / 3.0 * f[2]);

            const real_t alpha_0 = alpha_z(gamma_0, beta_0, tau_5, h[ida]);
            const real_t alpha_1 = alpha_z(gamma_1, beta_1, tau_5, h[ida]);
            const real_t alpha_2 = alpha_z(gamma_2, beta_2, tau_5, h[ida]);
            const real_t denom   = 1.0 / (alpha_0 + alpha_1 + alpha_2);
            const real_t w0      = alpha_0 * denom;
            const real_t w1      = alpha_1 * denom;
            const real_t w2      = alpha_2 * denom;

            const real_t fvel = 0.5 * (vel[m_idx_delta(-1, ida)] + vel[0]);            // i - 1/2
            flux[0]           = (m_sign(fvel) < 0.0) * (w0 * s0 + w1 * s1 + w2 * s2);  // it's a negative flux
        }
        // case 1: inteface i+1/2, left biased stencils - flux is positive
        {
            const real_t gamma_0 = 1.0 / 10.0;
            const real_t s0      = (1.0 / 3.0 * f[-2] - 7.0 / 6.0 * f[-1] + 11.0 / 6.0 * f[0]);
            const real_t gamma_1 = 3.0 / 5.0;
            const real_t s1      = (-1.0 / 6.0 * f[-1] + 5.0 / 6.0 * f[0] + 1.0 / 3.0 * f[1]);
            const real_t gamma_2 = 3.0 / 10.0;
            const real_t s2      = (1.0 / 3.0 * f[0] + 5.0 / 6.0 * f[1] - 1.0 / 6.0 * f[2]);

            const real_t alpha_0 = alpha_z(gamma_0, beta_0, tau_5, h[ida]);
            const real_t alpha_1 = alpha_z(gamma_1, beta_1, tau_5, h[ida]);
            const real_t alpha_2 = alpha_z(gamma_2, beta_2, tau_5, h[ida]);
            const real_t denom   = 1.0 / (alpha_0 + alpha_1 + alpha_2);
            const real_t w0      = alpha_0 * denom;
            const real_t w1      = alpha_1 * denom;
            const real_t w2      = alpha_2 * denom;

            const real_t fvel = 0.5 * (vel[m_idx_delta(+1, ida)] + vel[0]);            // i + 1/2
            flux[1]           = (m_sign(fvel) > 0.0) * (w0 * s0 + w1 * s1 + w2 * s2);  // it's a positive flux
        }
    };
    //-------------------------------------------------------------------------
    // apply the flux computations
    const real_t  alpha       = (accumulate_) ? 1.0 : 0.0;
    const real_t  oneoh[3]    = {1.0 / block->hgrid(0), 1.0 / block->hgrid(1), 1.0 / block->hgrid(2)};
    const real_t* data_vel[3] = {block->data(u_, 0).Read(), block->data(u_, 1).Read(), block->data(u_, 2).Read()};
    DoMagic_WENOZ(&flux_weno_5,
                  data_vel, block->data(fid_src, ida_).Read(), block->data(fid_trg, ida_).Write(),
                  is_outer, alpha, oneoh);
    // -------------------------------------------------------------------------
}
