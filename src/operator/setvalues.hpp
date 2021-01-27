#ifndef SRC_SETVALUES_HPP_
#define SRC_SETVALUES_HPP_

#include "core/doop.hpp"
#include "core/macros.hpp"
#include "core/pointers.hpp"
#include "core/types.hpp"
#include "grid/forestgrid.hpp"
#include "operator/blockoperator.hpp"
#include "prof.hpp"

//=====================================================================================================
/**
 * @brief defines a default SetValue object
 * 
 */
class SetValue : public BlockOperator {
   protected:
    lda_t  ida_start_ = 0;
    lda_t  ida_end_   = 0;
    real_t time_      = 0.0;

   public:
    explicit SetValue(m_ptr<const Wavelet> interp);

    void SetTime(const real_t time) { time_ = time; }

    void operator()(m_ptr<const ForestGrid> grid, m_ptr<Field> field);
    void operator()(m_ptr<const ForestGrid> grid, m_ptr<Field> field, const lda_t ida);

    /**
     * @brief Fill the value within the 
     * 
     * @param qid with be nullptr as it shouldn't be used here
     */
    virtual void FillGridBlock(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, m_ptr<Field> fid) = 0;
};

//=====================================================================================================
/**
 * @brief sets a the absolute value function in every dimension of a field
 */
class SetAbs : public SetValue {
   protected:
    real_t center_[3] = {0, 0, 0};
    real_t alpha_[3]  = {0.0, 0.0, 0.0};

    void FillGridBlock(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, m_ptr<Field> fid) override;

   public:
    SetAbs(const real_t alpha[3], const real_t center[3]);
    SetAbs(const real_t alpha[3], const real_t center[3], m_ptr<const Wavelet> interp);
};

//=====================================================================================================
/**
 * @brief sets a sinus in every dimension of a field:
 *
 *  formula: sum_i alpha[i] * sin(2*pi*freq[i] * x/L[i])
 * 
 */
class SetSinus : public SetValue {
   protected:
    real_t freq_[3]   = {0, 0, 0};
    real_t length_[3] = {0.0, 0.0, 0.0};
    real_t alpha_[3]  = {0.0, 0.0, 0.0};

    void FillGridBlock(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, m_ptr<Field> fid) override;

   public:
    SetSinus(const real_t length[3], const real_t freq[3], const real_t alpha[3]);
    SetSinus(const real_t length[3], const real_t freq[3], const real_t alpha[3], m_ptr<const Wavelet> interp);
};

//=====================================================================================================
/**
 * @brief sets a cosinus in every dimension of a field
 * 
 * formula: sum_i alpha[i] * cos(2*pi*freq[i] * x/L[i])
 * 
 */
class SetCosinus : public SetValue {
   protected:
    real_t freq_[3]   = {0, 0, 0};
    real_t length_[3] = {0.0, 0.0, 0.0};
    real_t alpha_[3]  = {0.0, 0.0, 0.0};

    void FillGridBlock(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, m_ptr<Field> fid) override;

   public:
    SetCosinus(const real_t length[3], const real_t freq[3], const real_t alpha[3]);
    SetCosinus(const real_t length[3], const real_t freq[3], const real_t alpha[3], m_ptr<const Wavelet> interp);
};

//=====================================================================================================
/**
 * @brief sets a polynomial in every dimension of a field
 */
class SetPolynom : public SetValue {
   protected:
    bool   extend_   = false;            //!< indicate if we need to fill the ghosts or not
    lid_t  deg_[3]   = {0, 0, 0};        //!< the degree of the polynomial
    real_t dir_[3]   = {0.0, 0.0, 0.0};  //!< the direction concerned: 1.0 means involved, 0.0 means not involved
    real_t shift_[3] = {0.0, 0.0, 0.0};

    void FillGridBlock(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, m_ptr<Field> fid) override;

   public:
    SetPolynom(const lid_t degree[3], const real_t direction[3], const real_t shift[3]);
    SetPolynom(const lid_t degree[3], const real_t direction[3], const real_t shift[3], m_ptr<const Wavelet> interp);
};

//=====================================================================================================
/**
 * @brief set an exponential, centered around center_ and with a sigma = sigma_ and integral = alpha_
 */
class SetExponential : public SetValue {
   protected:
    real_t center_[3] = {0.0, 0.0, 0.0};
    real_t sigma_[3]  = {0.0, 0.0, 0.0};
    real_t alpha_     = 1.0;

    void FillGridBlock(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, m_ptr<Field> fid) override;

   public:
    SetExponential(const real_t center[3], const real_t sigma[3], const real_t alpha);
    SetExponential(const real_t center[3], const real_t sigma[3], const real_t alpha, m_ptr<const Wavelet> interp);
};

//=====================================================================================================
class SetErf : public SetValue {
   protected:
    real_t center_[3] = {0.0, 0.0, 0.0};
    real_t sigma_[3]  = {0.0, 0.0, 0.0};
    real_t alpha_     = 1.0;

    void FillGridBlock(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, m_ptr<Field> fid) override;

   public:
    SetErf(const real_t center[3], const real_t sigma[3], const real_t alpha);
    SetErf(const real_t center[3], const real_t sigma[3], const real_t alpha, m_ptr<const Wavelet> interp);
};

//=====================================================================================================
class SetVortexRing : public SetValue {
   protected:
    lda_t  normal_    = 0;                //!< the direction normal to the ring, i.e. the z direction
    real_t sigma_     = 0.0;              //!< the direction normal to the ring, i.e. the z direction
    real_t radius_    = 0.0;              //!< the direction normal to the ring, i.e. the z direction
    real_t center_[3] = {0.0, 0.0, 0.0};  //!< the center of the ring

    void FillGridBlock(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, m_ptr<Field> fid) override;

   public:
    SetVortexRing(const lda_t normal, const real_t center[3], const real_t sigma, const real_t radius);
    SetVortexRing(const lda_t normal, const real_t center[3], const real_t sigma, const real_t radius, m_ptr<const Wavelet> interp);
};

//=====================================================================================================
class SetCompactVortexRing : public SetValue {
   protected:
    lda_t  normal_    = 0;                //!< the direction normal to the ring, i.e. the z direction
    real_t sigma_     = 0.0;              //!< the direction normal to the ring, i.e. the z direction
    real_t radius_    = 0.0;              //!< the direction normal to the ring, i.e. the z direction
    real_t center_[3] = {0.0, 0.0, 0.0};  //!< the center of the ring
    real_t cutoff_    = 0.0;              //!< the cutoff distance, i.e. the distance after which the gaussian is set to 0.0

    void FillGridBlock(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, m_ptr<Field> fid) override;

   public:
    SetCompactVortexRing(const lda_t normal, const real_t center[3], const real_t sigma, const real_t radius, const real_t cutoff);
    SetCompactVortexRing(const lda_t normal, const real_t center[3], const real_t sigma, const real_t radius, const real_t cutoff, m_ptr<const Wavelet> interp);
};

//=====================================================================================================
class SetABSVelocity : public SetValue {
   protected:
    real_t a_ = 0;    //!< the direction normal to the ring, i.e. the z direction
    real_t b_ = 0.0;  //!< the direction normal to the ring, i.e. the z direction
    real_t c_ = 0.0;  //!< the direction normal to the ring, i.e. the z direction

    void FillGridBlock(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, m_ptr<Field> fid) override;

   public:
    SetABSVelocity(const real_t a, const real_t b, const real_t c);
    SetABSVelocity(const real_t a, const real_t b, const real_t c, m_ptr<const Wavelet> interp);
};

//=====================================================================================================
class SetScalarRing : public SetValue {
   protected:
    lda_t  normal_    = 0;                //!< the direction normal to the ring, i.e. the z direction
    real_t sigma_     = 0.0;              //!< the direction normal to the ring, i.e. the z direction
    real_t radius_    = 0.0;              //!< the direction normal to the ring, i.e. the z direction
    real_t center_[3] = {0.0, 0.0, 0.0};  //!< the center of the ring
    real_t vel_[3]    = {0.0, 0.0, 0.0};  //!< advection velocity for the ring

    void FillGridBlock(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, m_ptr<Field> fid) override;

   public:
    SetScalarRing(const lda_t normal, const real_t center[3], const real_t sigma, const real_t radius, const real_t vel[3]);
    SetScalarRing(const lda_t normal, const real_t center[3], const real_t sigma, const real_t radius, const real_t vel[3], m_ptr<const Wavelet> interp);
};

//=====================================================================================================
class SetScalarTube : public SetValue {
   protected:
    lda_t  dir_       = 0;                //!< the direction normal of the tube
    real_t sigma_     = 0.0;              //!< the sigma of one tube
    real_t b_         = 0.0;              //!< the distance between two tubes
    real_t center_[3] = {0.0, 0.0, 0.0};  //!< the center of the tube system

    void FillGridBlock(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, m_ptr<Field> fid) override;

   public:
    SetScalarTube(const lda_t dir, const real_t center[3], const real_t sigma, const real_t b);
    SetScalarTube(const lda_t dir, const real_t center[3], const real_t sigma, const real_t b, m_ptr<const Wavelet> interp);
};

#endif  // SRC_SETVALUES_HPP_
