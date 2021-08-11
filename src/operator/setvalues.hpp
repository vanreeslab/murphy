#ifndef SRC_SETVALUES_HPP_
#define SRC_SETVALUES_HPP_

#include <vector>

#include "core/doop.hpp"
#include "core/forloop.hpp"
#include "core/macros.hpp"
#include "core/types.hpp"
#include "grid/field.hpp"
#include "grid/forestgrid.hpp"
#include "operator/blockoperator.hpp"
#include "tools/prof.hpp"

// using lambda_setvalue_scalar_pos_t = lambda_t<void, const real_t*, real_t*>;
// using lambda_setvalue_vector_pos_t = lambda_t<void, const real_t*, real_t**>;
using lambda_setvalue_t = lambda_i3_t<void, const CartBlock* const, const Field* const>;

//=====================================================================================================
/**
 * @brief Applies the lambda_i3block_t expr to the block
 * 
 */
class SetValue : public BlockOperator {
   private:
    const lambda_setvalue_t& expr_;

   public:
    explicit SetValue() = delete;
    explicit SetValue(const lambda_setvalue_t& expr, const bidx_t* ghost_len = nullptr) : expr_(expr), BlockOperator(ghost_len){};

    /**
     * @brief set the expr_ to the grid in the given field
     * 
     * @warning the field could have been captured by the lambda but it's less general and prevents the assignation of the ghost status
     * 
     * @param grid 
     * @param field 
     */
    void operator()(const ForestGrid* grid, Field* const field) const {
        m_begin;
        //-------------------------------------------------------------------------
        // go for it
        DoOpMesh(this, &SetValue::FillGridBlock, grid, field);
        // update the ghost status
        field->ghost_len(ghost_len_res_);
        //-------------------------------------------------------------------------
        m_end;
    };
    void FillGridBlock(const qid_t* qid, CartBlock* block, const Field* const fid) const {
        //-------------------------------------------------------------------------
        // translate the expr into a know pattern for the forloop
        auto op = [=](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
            this->expr_(i0, i1, i2, block, fid);
        };
        // run the forloop on it
        for_loop(&op, start_, end_);
        //-------------------------------------------------------------------------
    };
};

// template <>
// inline void SetValue<lambda_setvalue_idx_t>::FillGridBlock(const qid_t* qid, CartBlock* block, const Field* const fid) const {
//     //-------------------------------------------------------------------------
//
//     //-------------------------------------------------------------------------
// };
// template <>
// inline void SetValue<lambda_setvalue_scalar_pos_t>::FillGridBlock(const qid_t* qid, CartBlock* block, const Field* const fid) const {
//     //-------------------------------------------------------------------------
//     // translate the expr into a know pattern for the forloop
//     auto op = [=](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
//         // get the position
//         real_t pos[3];
//         block->pos(i0, i1, i2, pos);

//         real_t* data = block->data(fid).Write(i0, i1, i2);
//         expr_(pos, data);
//     };
//     // run the forloop on it
//     for_loop(&op, start_, end_);
//     //-------------------------------------------------------------------------
// };
// template <>
// inline void SetValue<lambda_setvalue_vector_pos_t>::FillGridBlock(const qid_t* qid, CartBlock* block, const Field* const fid) const {

// };

// declare later defined functions as lambdas
extern lambda_t<real_t, const real_t[], const real_t[], const real_t>                            scalar_exp;
extern lambda_t<real_t, const real_t[], const real_t[], const real_t, const lda_t>               scalar_tube;
extern lambda_t<real_t, const real_t[], const real_t[], const real_t, const real_t, const lda_t> scalar_ring;
extern lambda_t<real_t, const real_t[], const real_t[],
                const lda_t, const real_t, const real_t, const real_t, const std::vector<short_t>, const std::vector<real_t> >
    scalar_compact_ring;

//  = [](const bidx_t i0, const bidx_t i1, const bidx_t i2, const CartBlock* const block) -> real_t {
//     // get the position
//     real_t pos[3];
//     m_pos(pos, i0, i1, i2, block->hgrid(), block->xyz());

//     const real_t sigma     = 0.05;
//     const real_t center[3] = {0.5, 0.5, 0.5};

//     // compute the gaussian
//     const real_t rhox = (pos[0] - center[0]) / sigma;
//     const real_t rhoy = (pos[1] - center[1]) / sigma;
//     const real_t rhoz = (pos[2] - center[2]) / sigma;
//     const real_t rho  = rhox * rhox + rhoy * rhoy + rhoz * rhoz;

//     return std::exp(-rho);
// };

// //=====================================================================================================
// /**
//  * @brief sets a the absolute value function in every dimension of a field
//  */
// class SetAbs : public SetValue {
//    protected:
//     real_t center_[3] = {0, 0, 0};
//     real_t alpha_[3]  = {0.0, 0.0, 0.0};

//     void FillGridBlock(const qid_t*  qid, GridBlock*  block, Field*  fid) override;

//    public:
//     SetAbs(const real_t alpha[3], const real_t center[3]);
//     SetAbs(const real_t alpha[3], const real_t center[3], const Wavelet*  interp);
// };

// //=====================================================================================================
// /**
//  * @brief sets a sinus in every dimension of a field:
//  *
//  *  formula: sum_i alpha[i] * sin(2*pi*freq[i] * x/L[i])
//  *
//  */
// class SetSinus : public SetValue {
//    protected:
//     real_t freq_[3]   = {0, 0, 0};
//     real_t length_[3] = {0.0, 0.0, 0.0};
//     real_t alpha_[3]  = {0.0, 0.0, 0.0};

//     void FillGridBlock(const qid_t*  qid, GridBlock*  block, Field*  fid) override;

//    public:
//     SetSinus(const real_t length[3], const real_t freq[3], const real_t alpha[3]);
//     SetSinus(const real_t length[3], const real_t freq[3], const real_t alpha[3], const Wavelet*  interp);
// };

// //=====================================================================================================
// /**
//  * @brief sets a cosinus in every dimension of a field
//  *
//  * formula: sum_i alpha[i] * cos(2*pi*freq[i] * x/L[i])
//  *
//  */
// class SetCosinus : public SetValue {
//    protected:
//     real_t freq_[3]   = {0, 0, 0};
//     real_t length_[3] = {0.0, 0.0, 0.0};
//     real_t alpha_[3]  = {0.0, 0.0, 0.0};

//     void FillGridBlock(const qid_t*  qid, GridBlock*  block, Field*  fid) override;

//    public:
//     SetCosinus(const real_t length[3], const real_t freq[3], const real_t alpha[3]);
//     SetCosinus(const real_t length[3], const real_t freq[3], const real_t alpha[3], const Wavelet*  interp);
// };

// //=====================================================================================================
// /**
//  * @brief sets a polynomial in every dimension of a field
//  */
// class SetPolynom : public SetValue {
//    protected:
//     bool   extend_   = false;            //!< indicate if we need to fill the ghosts or not
//     lid_t  deg_[3]   = {0, 0, 0};        //!< the degree of the polynomial
//     real_t dir_[3]   = {0.0, 0.0, 0.0};  //!< the direction concerned: 1.0 means involved, 0.0 means not involved
//     real_t shift_[3] = {0.0, 0.0, 0.0};

//     void FillGridBlock(const qid_t*  qid, GridBlock*  block, Field*  fid) override;

//    public:
//     SetPolynom(const lid_t degree[3], const real_t direction[3], const real_t shift[3]);
//     SetPolynom(const lid_t degree[3], const real_t direction[3], const real_t shift[3], const Wavelet*  interp);
// };

// //=====================================================================================================
// /**
//  * @brief set an exponential, centered around center_ and with a sigma = sigma_ and integral = alpha_
//  */
// class SetExponential : public SetValue {
//    protected:
//     real_t center_[3] = {0.0, 0.0, 0.0};
//     real_t sigma_[3]  = {0.0, 0.0, 0.0};
//     real_t alpha_     = 1.0;

//     void FillGridBlock(const qid_t*  qid, GridBlock*  block, Field*  fid) override;

//    public:
//     SetExponential(const real_t center[3], const real_t sigma[3], const real_t alpha);
//     SetExponential(const real_t center[3], const real_t sigma[3], const real_t alpha, const Wavelet*  interp);
// };

// //=====================================================================================================
// class SetErf : public SetValue {
//    protected:
//     real_t center_[3] = {0.0, 0.0, 0.0};
//     real_t sigma_[3]  = {0.0, 0.0, 0.0};
//     real_t alpha_     = 1.0;

//     void FillGridBlock(const qid_t*  qid, GridBlock*  block, Field*  fid) override;

//    public:
//     SetErf(const real_t center[3], const real_t sigma[3], const real_t alpha);
//     SetErf(const real_t center[3], const real_t sigma[3], const real_t alpha, const Wavelet*  interp);
// };

// //=====================================================================================================
// class SetVortexRing : public SetValue {
//    protected:
//     lda_t  normal_    = 0;                //!< the direction normal to the ring, i.e. the z direction
//     real_t sigma_     = 0.0;              //!< the direction normal to the ring, i.e. the z direction
//     real_t radius_    = 0.0;              //!< the direction normal to the ring, i.e. the z direction
//     real_t center_[3] = {0.0, 0.0, 0.0};  //!< the center of the ring

//     void FillGridBlock(const qid_t*  qid, GridBlock*  block, Field*  fid) override;

//    public:
//     SetVortexRing(const lda_t normal, const real_t center[3], const real_t sigma, const real_t radius);
//     SetVortexRing(const lda_t normal, const real_t center[3], const real_t sigma, const real_t radius, const Wavelet*  interp);
// };

// //=====================================================================================================
// class SetCompactVortexRing : public SetValue {
//    protected:
//     lda_t  normal_    = 0;                //!< the direction normal to the ring, i.e. the z direction
//     real_t sigma_     = 0.0;              //!< the direction normal to the ring, i.e. the z direction
//     real_t radius_    = 0.0;              //!< the direction normal to the ring, i.e. the z direction
//     real_t center_[3] = {0.0, 0.0, 0.0};  //!< the center of the ring
//     real_t cutoff_    = 0.0;              //!< the cutoff distance, i.e. the distance after which the gaussian is set to 0.0

//     void FillGridBlock(const qid_t*  qid, GridBlock*  block, Field*  fid) override;

//    public:
//     SetCompactVortexRing(const lda_t normal, const real_t center[3], const real_t sigma, const real_t radius, const real_t cutoff);
//     SetCompactVortexRing(const lda_t normal, const real_t center[3], const real_t sigma, const real_t radius, const real_t cutoff, const Wavelet*  interp);
// };

// //=====================================================================================================
// class SetABSVelocity : public SetValue {
//    protected:
//     real_t a_ = 0;    //!< the direction normal to the ring, i.e. the z direction
//     real_t b_ = 0.0;  //!< the direction normal to the ring, i.e. the z direction
//     real_t c_ = 0.0;  //!< the direction normal to the ring, i.e. the z direction

//     void FillGridBlock(const qid_t*  qid, GridBlock*  block, Field*  fid) override;

//    public:
//     SetABSVelocity(const real_t a, const real_t b, const real_t c);
//     SetABSVelocity(const real_t a, const real_t b, const real_t c, const Wavelet*  interp);
// };

// //=====================================================================================================
// class SetScalarRing : public SetValue {
//    protected:
//     lda_t  normal_    = 0;                //!< the direction normal to the ring, i.e. the z direction
//     real_t sigma_     = 0.0;              //!< the direction normal to the ring, i.e. the z direction
//     real_t radius_    = 0.0;              //!< the direction normal to the ring, i.e. the z direction
//     real_t center_[3] = {0.0, 0.0, 0.0};  //!< the center of the ring
//     real_t vel_[3]    = {0.0, 0.0, 0.0};  //!< advection velocity for the ring

//     void FillGridBlock(const qid_t*  qid, GridBlock*  block, Field*  fid) override;

//    public:
//     SetScalarRing(const lda_t normal, const real_t center[3], const real_t sigma, const real_t radius, const real_t vel[3]);
//     SetScalarRing(const lda_t normal, const real_t center[3], const real_t sigma, const real_t radius, const real_t vel[3], const Wavelet*  interp);
// };

// //=====================================================================================================
// class SetScalarTube : public SetValue {
//    protected:
//     lda_t  dir_       = 0;                //!< the direction normal of the tube
//     real_t sigma_     = 0.0;              //!< the sigma of one tube
//     real_t b_         = 0.0;              //!< the distance between two tubes
//     real_t center_[3] = {0.0, 0.0, 0.0};  //!< the center of the tube system

//     void FillGridBlock(const qid_t*  qid, GridBlock*  block, Field*  fid) override;

//    public:
//     SetScalarTube(const lda_t dir, const real_t center[3], const real_t sigma, const real_t b);
//     SetScalarTube(const lda_t dir, const real_t center[3], const real_t sigma, const real_t b, const Wavelet*  interp);
// };

#endif  // SRC_SETVALUES_HPP_
