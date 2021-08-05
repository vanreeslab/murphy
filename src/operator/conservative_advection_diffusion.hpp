// #ifndef SRC_CONSERVATIVE_ADVECTION_DIFFUSION_HPP_
// #define SRC_CONSERVATIVE_ADVECTION_DIFFUSION_HPP_

// #include "core/macros.hpp"
// #include "core/types.hpp"
// #include "forloop.hpp"
// #include "stencil.hpp"

// /**
//  * @brief compute the advection + diffusion stencil on a given field
//  * 
//  * @tparam length_advection the size of the advection stencil. (e.g. 3 for 2nd order and 5 for 4th order)
//  * @tparam length_diffusion the size of the diffusion stencil. (e.g. 3 for 2nd order and 5 for 4th order)
//  */
// template <sid_t length_advection, sid_t length_diffusion>
// class Conservative_AdvectionDiffusion : public Stencil {
//    protected:
//     real_t nu_          = 0.0;
//     real_t u_stream_[3] = {0.0, 0.0, 1.0};

//    public:
//     explicit Conservative_AdvectionDiffusion(const real_t nu, const real_t u_stream[3]) : Stencil() {
//         //-------------------------------------------------------------------------
//         nu_ = nu;
//         for (lda_t i = 0; i < 3; ++i) {
//             u_stream_[i] = u_stream[i];
//         }
//         m_log("Conservative Advection Diffusion with u stream = %e %e %e and nu = %e", u_stream_[0], u_stream_[1], u_stream_[2], nu_);
//         //-------------------------------------------------------------------------
//     }

//     virtual lid_t NGhost() const override { return m_max(length_advection / 2, length_diffusion / 2); };

//    protected:
//     void DoMagic(const qid_t*  qid, GridBlock*  block, const bool is_outer, const Field*  fid_src, Field*  fid_trg) const override;
// };

// template <sid_t length_advection, sid_t length_diffusion>
// void Conservative_AdvectionDiffusion<length_advection, length_diffusion>::DoMagic(const qid_t*  qid, GridBlock*  block, const bool is_outer, const Field*  fid_src, Field*  fid_trg) {
//     m_assert((fid_src->lda() == 3) && (fid_trg->lda() == 3), "the source must be 3 times smaller than the target ");
//     static_assert(length_advection == 4, "the size of the stencil is not supported");
//     static_assert(length_diffusion == 3 || length_diffusion == 5, "the size of the stencil is not supported");
//     //-------------------------------------------------------------------------
//     const real_t scale_d[3]  = {1.0 / block->hgrid(0), 1.0 / block->hgrid(1), 1.0 / block->hgrid(2)};
//     const real_t scale_d2[3] = {nu_ / pow(block->hgrid(0), 2), nu_ / pow(block->hgrid(1), 2), nu_ / pow(block->hgrid(2), 2)};

//     const real_t* data_src = block->data(fid_src, ida_).Read();
//     real_t*       data_trg = block->data(fid_trg, ida_).Write();

//     constexpr real_t two_third  = 2.0 / 3.0;
//     constexpr real_t one_twelve = 1.0 / 12.0;
//     constexpr real_t four_third = 4.0 / 3.0;
//     constexpr real_t five_half  = 5.0 / 2.0;
//     constexpr real_t one_sixth  = 1.0 / 6.0;
//     constexpr real_t five_sixth = 5.0 / 6.0;
//     constexpr real_t one_third  = 1.0 / 3.0;

//     auto op = [=, &data_trg, *this](const lid_t i0, const lid_t i1, const lid_t i2) -> void {
//         // get the data pointer in front of the row for every cache line
//         real_t*       trg  = data_trg + m_idx(i0, i1, i2);                // cache line for writting
//         const real_t* src  = data_src + m_idx(i0, i1, i2);                // cache line for reading
//         const real_t  u[3] = {u_stream_[0], u_stream_[1], u_stream_[2]};  // todo add the real velocity here and do a proper u(+1/2) interpolation

//         //advection
//         // m_assert(trg[0] == 0.0,"the target value must be 0.0");
//         if constexpr (length_advection == 4) {
//             // X -
//             const real_t f_mx_p = (-one_sixth * src[m_idx(-2, 0, 0)] + five_sixth * src[m_idx(-1, 0, 0)] + one_third * src[m_idx(0, 0, 0)]);
//             const real_t f_mx_m = (one_third * src[m_idx(-1, 0, 0)] + five_sixth * src[m_idx(0, 0, 0)] - one_sixth * src[m_idx(+1, 0, 0)]);
//             // X +
//             const real_t f_px_p = (-one_sixth * src[m_idx(-1, 0, 0)] + five_sixth * src[m_idx(0, 0, 0)] + one_third * src[m_idx(+1, 0, 0)]);
//             const real_t f_px_m = (one_third * src[m_idx(0, 0, 0)] + five_sixth * src[m_idx(+1, 0, 0)] - one_sixth * src[m_idx(+2, 0, 0)]);
//             // Y -
//             const real_t f_my_p = (-one_sixth * src[m_idx(0, -2, 0)] + five_sixth * src[m_idx(0, -1, 0)] + one_third * src[m_idx(0, 0, 0)]);
//             const real_t f_my_m = (one_third * src[m_idx(0, -1, 0)] + five_sixth * src[m_idx(0, 0, 0)] - one_sixth * src[m_idx(0, +1, 0)]);
//             // Y +
//             const real_t f_py_p = (-one_sixth * src[m_idx(0, -1, 0)] + five_sixth * src[m_idx(0, 0, 0)] + one_third * src[m_idx(0, +1, 0)]);
//             const real_t f_py_m = (one_third * src[m_idx(0, 0, 0)] + five_sixth * src[m_idx(0, +1, 0)] - one_sixth * src[m_idx(0, +2, 0)]);
//             // Z -
//             const real_t f_mz_p = (-one_sixth * src[m_idx(0, 0, -2)] + five_sixth * src[m_idx(0, 0, -1)] + one_third * src[m_idx(0, 0, 0)]);
//             const real_t f_mz_m = (one_third * src[m_idx(0, 0, -1)] + five_sixth * src[m_idx(0, 0, 0)] - one_sixth * src[m_idx(0, 0, +1)]);
//             // Z +
//             const real_t f_pz_p = (-one_sixth * src[m_idx(0, 0, -1)] + five_sixth * src[m_idx(0, 0, 0)] + one_third * src[m_idx(0, 0, +1)]);
//             const real_t f_pz_m = (one_third * src[m_idx(0, 0, 0)] + five_sixth * src[m_idx(0, 0, +1)] - one_sixth * src[m_idx(0, 0, +2)]);

//             // add
//             (*trg) -= scale_d[0] * (m_max(0.0, u[0]) * (f_px_p - f_mx_p) + m_min(0.0, u[0]) * (f_px_m - f_mx_m));
//             (*trg) -= scale_d[1] * (m_max(0.0, u[1]) * (f_py_p - f_my_p) + m_min(0.0, u[1]) * (f_py_m - f_my_m));
//             (*trg) -= scale_d[2] * (m_max(0.0, u[2]) * (f_pz_p - f_mz_p) + m_min(0.0, u[2]) * (f_pz_m - f_mz_m));
//         }
//         // if (pos[0] == 0.75 && pos[1] == 0.75 && pos[2] == 0.5 && i0 == 1 && i1 == 1 && i2 == 0) {
//         //     // m_log("block num %d", qid->cid);
//         //     m_log("block num %d (ida = %d): trg = %e, src = %e", qid->cid, ida_, (*trg), (*src));
//         // }
//         // if (pos[0] == 0.125 && pos[1] == 0.125 && pos[2] == 0.5 && i0 == 15 && i1 == 15 && i2 == 0) {
//         //     m_log("block num %d (ida = %d): trg = %e, src = %e", qid->cid, ida_, (*trg), (*src));
//         // }
//         // if (pos[0] == 0.75 && pos[1] == 0.75 && pos[2] == 0.5 && i0 == 0 && i1 == 0 && i2 == 0) {
//         //     // m_log("block num %d", qid->cid);
//         //     m_log("block num %d (ida = %d): trg = %e, src = %e", qid->cid, ida_, (*trg), (*src));
//         // }
//         // if (pos[0] == 0.25 && pos[1] == 0.25 && pos[2] == 0.5 && i0 == 0 && i1 == 0 && i2 == 0) {
//         //     m_log("block num %d (ida = %d): trg = %e, src = %e", qid->cid, ida_, (*trg), (*src));
//         // }
//         // diffusion
//         if constexpr (length_diffusion == 3) {
//             (*trg) += scale_d2[0] * (src[m_idx(+1, 0, 0)] - 2.0 * src[m_idx(0, 0, 0)] + src[m_idx(-1, 0, 0)]);
//             (*trg) += scale_d2[1] * (src[m_idx(0, +1, 0)] - 2.0 * src[m_idx(0, 0, 0)] + src[m_idx(0, -1, 0)]);
//             (*trg) += scale_d2[2] * (src[m_idx(0, 0, +1)] - 2.0 * src[m_idx(0, 0, 0)] + src[m_idx(0, 0, -1)]);
//         } else if constexpr (length_diffusion == 5) {
//             (*trg) += scale_d2[0] * (-one_twelve * src[m_idx(+2, 0, 0)] + four_third * src[m_idx(+1, 0, 0)] - five_half * src[m_idx(0, 0, 0)] + four_third * src[m_idx(-1, 0, 0)] - one_twelve * src[m_idx(-2, 0, 0)]);
//             (*trg) += scale_d2[1] * (-one_twelve * src[m_idx(0, +2, 0)] + four_third * src[m_idx(0, +1, 0)] - five_half * src[m_idx(0, 0, 0)] + four_third * src[m_idx(0, -1, 0)] - one_twelve * src[m_idx(0, -2, 0)]);
//             (*trg) += scale_d2[2] * (-one_twelve * src[m_idx(0, 0, +2)] + four_third * src[m_idx(0, 0, +1)] - five_half * src[m_idx(0, 0, 0)] + four_third * src[m_idx(0, 0, -1)] - one_twelve * src[m_idx(0, 0, -2)]);
//         }
//     };

//     if (!is_outer) {
//         for_loop<M_GS, M_N - M_GS>(&op);
//     } else {
//         // do the most on the X side to use vectorization
//         for_loop<0, M_N, 0, M_N, 0, M_GS>(&op);          // Z-
//         for_loop<0, M_N, 0, M_N, M_N - M_GS, M_N>(&op);  // Z+

//         for_loop<0, M_N, 0, M_GS, M_GS, M_N - M_GS>(&op);          // Y-
//         for_loop<0, M_N, M_N - M_GS, M_N, M_GS, M_N - M_GS>(&op);  // Y+

//         for_loop<0, M_GS, M_GS, M_N - M_GS, M_GS, M_N - M_GS>(&op);          // X-
//         for_loop<M_N - M_GS, M_N, M_GS, M_N - M_GS, M_GS, M_N - M_GS>(&op);  // X+
//     }
//     //-------------------------------------------------------------------------
// }

// #endif  // SRC_ADVECTION_DIFFUSION_HPP_