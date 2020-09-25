#ifndef SRC_SETVALUES_HPP_
#define SRC_SETVALUES_HPP_

#include "doop.hpp"
#include "forestgrid.hpp"
#include "murphy.hpp"

//=====================================================================================================
/**
 * @brief defines a default SetValue object
 * 
 */
class SetValue {
   protected:
    lid_t start_ = 0;    //!< the starting index in 3D = [start_,start_,start_]
    lid_t end_   = M_N;  //!< the ending index in 3D = [end_,end_,end_]

   public:
    explicit SetValue(const lid_t nghost_front, const lid_t nghost_back);

    virtual void operator()(ForestGrid* grid, Field* field);
    virtual void FillGridBlock(const qid_t* qid, GridBlock* block, Field* fid) = 0;
};

//=====================================================================================================
/**
 * @brief sets a the absolute value function in every dimension of a field
 */
class SetAbs : public SetValue {
   protected:
    real_t center_[3] = {0, 0, 0};
    real_t alpha_[3]  = {0.0, 0.0, 0.0};

    void FillGridBlock(const qid_t* qid, GridBlock* block, Field* fid) override;

   public:
    SetAbs(const real_t alpha[3], const real_t center[3]);
    SetAbs(const real_t alpha[3], const real_t center[3], const lid_t nghost_front, const lid_t nghost_back);
};

//=====================================================================================================
/**
 * @brief sets a sinus in every dimension of a field
 */
class SetSinus : public SetValue {
   protected:
    real_t freq_[3]   = {0, 0, 0};
    real_t length_[3] = {0.0, 0.0, 0.0};

    void FillGridBlock(const qid_t* qid, GridBlock* block, Field* fid) override;

   public:
    SetSinus(const real_t length[3], const real_t freq[3]);
    SetSinus(const real_t length[3], const real_t freq[3], const lid_t nghost_front, const lid_t nghost_back);
};

//=====================================================================================================
/**
 * @brief sets a cosinus in every dimension of a field
 */
class SetCosinus : public SetValue {
   protected:
    real_t freq_[3]   = {0, 0, 0};
    real_t length_[3] = {0.0, 0.0, 0.0};

    void FillGridBlock(const qid_t* qid, GridBlock* block, Field* fid) override;

   public:
    SetCosinus(const real_t length[3], const real_t freq[3]);
    SetCosinus(const real_t length[3], const real_t freq[3], const lid_t nghost_front, const lid_t nghost_back);
};

//=====================================================================================================
/**
 * @brief sets a polynomial in every dimension of a field
 */
class SetPolynom : public SetValue {
   protected:
    bool   extend_ = false;            //!< indicate if we need to fill the ghosts or not
    lid_t  deg_[3] = {0, 0, 0};        //!< the degree of the polynomial
    real_t dir_[3] = {0.0, 0.0, 0.0};  //!< the direction concerned: 1.0 means involved, 0.0 means not involved

    void FillGridBlock(const qid_t* qid, GridBlock* block, Field* fid) override;

   public:
    SetPolynom(const lid_t degree[3], const real_t direction[3]);
    SetPolynom(const lid_t degree[3], const real_t direction[3], const lid_t nghost_front, const lid_t nghost_back);
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

    void FillGridBlock(const qid_t* qid, GridBlock* block, Field* fid) override;

   public:
    SetExponential(const real_t center[3], const real_t sigma[3], const real_t alpha);
    SetExponential(const real_t center[3], const real_t sigma[3], const real_t alpha, const lid_t nghost_front, const lid_t nghost_back);
};

//=====================================================================================================
class SetErf : public SetValue {
   protected:
    real_t center_[3] = {0.0, 0.0, 0.0};
    real_t sigma_[3]  = {0.0, 0.0, 0.0};
    real_t alpha_     = 1.0;

    void FillGridBlock(const qid_t* qid, GridBlock* block, Field* fid) override;

   public:
    SetErf(const real_t center[3], const real_t sigma[3], const real_t alpha);
    SetErf(const real_t center[3], const real_t sigma[3], const real_t alpha, const lid_t nghost_front, const lid_t nghost_back);
};

//=====================================================================================================
class SetVortexRing : public SetValue {
   protected:
    lda_t  normal_    = 0;                //!< the direction normal to the ring, i.e. the z direction
    real_t sigma_     = 0.0;              //!< the direction normal to the ring, i.e. the z direction
    real_t radius_    = 0.0;              //!< the direction normal to the ring, i.e. the z direction
    real_t center_[3] = {0.0, 0.0, 0.0};  //!< the center of the ring

    void FillGridBlock(const qid_t* qid, GridBlock* block, Field* fid) override;

   public:
    SetVortexRing(const lda_t normal, const real_t center[3], const real_t sigma, const real_t radius);
    SetVortexRing(const lda_t normal, const real_t center[3], const real_t sigma, const real_t radius, const lid_t nghost_front, const lid_t nghost_back);
};

#endif  // SRC_SETVALUES_HPP_
