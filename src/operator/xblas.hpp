#ifndef SRC_OPERATOR_XBLAS_HPP
#define SRC_OPERATOR_XBLAS_HPP


#include "core/pointers.hpp"
#include "operator/blockoperator.hpp"

/**
 * @brief perform the max(fabs()) operation on a block, i.e. return the infinite norm of a field
 *
 */
class BMax : public BlockOperator {
   protected:
    real_t max_;
    
   public:
    explicit BMax() noexcept;
    explicit BMax(const Wavelet*  interp) noexcept;

    real_t operator()(const ForestGrid*  grid, const Field*  fid_x);
    void   ComputeBMaxGridBlock(const qid_t*  qid, GridBlock*  block, const Field*  fid_x);
};

/**
 * @brief perform the min and max operation on a block and return the two values after and AllReduce
 *
 */
class BMinMax : public BlockOperator {
   protected:
    lda_t ida_;
    real_t max_, min_;

   public:
    explicit BMinMax() noexcept;
    explicit BMinMax(const Wavelet*  interp) noexcept;

    void operator()(const ForestGrid*  grid, const Field*  fid_x, real_t* min, real_t* max);
    void ComputeBMinMaxGridBlock(const qid_t*  qid, GridBlock*  block, const Field*  fid_x);
};

class BMoment : public BlockOperator {
   protected:
    lda_t  ida_;
    real_t moment0_;
    real_t moment1_[3];

   public:
    explicit BMoment() noexcept;
    explicit BMoment(const Wavelet*  interp) noexcept;

    void operator()(const ForestGrid*  grid, const Field*  fid_x, real_t* moment0, real_t* moment1);
    void ComputeBMomentGridBlock(const qid_t*  qid, GridBlock*  block, const Field*  fid_x);
};

class BMean : public BlockOperator {
   protected:
    lda_t  ida_;
    real_t sum_;

   public:
    explicit BMean() noexcept;
    explicit BMean(const Wavelet*  interp) noexcept;

    void operator()(const ForestGrid*  grid, const Field*  fid_x, real_t* sum);
    void ComputeBMeanGridBlock(const qid_t*  qid, GridBlock*  block, const Field*  fid_x);
};

class BDiscreteMoment : public BlockOperator {
   protected:
    lda_t  ida_;
    real_t moment0_;
    real_t moment1_[3];

   public:
    explicit BDiscreteMoment() noexcept;
    explicit BDiscreteMoment(const Wavelet*  interp) noexcept;

    void operator()(const ForestGrid*  grid, const Field*  fid_x, real_t* moment0, real_t* moment1);
    void ComputeBDiscreteMomentGridBlock(const qid_t*  qid, GridBlock*  block, const Field*  fid_x);
};

#endif