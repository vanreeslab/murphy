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
    explicit BMax(m_ptr<const Wavelet> interp) noexcept;

    real_t operator()(m_ptr<const ForestGrid> grid, m_ptr<const Field> fid_x);
    void   ComputeBMaxGridBlock(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, m_ptr<const Field> fid_x);
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
    explicit BMinMax(m_ptr<const Wavelet> interp) noexcept;

    void operator()(m_ptr<const ForestGrid> grid, m_ptr<const Field> fid_x, real_t* min, real_t* max);
    void ComputeBMinMaxGridBlock(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, m_ptr<const Field> fid_x);
};

class BMoment : public BlockOperator {
   protected:
    lda_t  ida_;
    real_t moment0_;
    real_t moment1_[3];

   public:
    explicit BMoment() noexcept;
    explicit BMoment(m_ptr<const Wavelet> interp) noexcept;

    void operator()(m_ptr<const ForestGrid> grid, m_ptr<const Field> fid_x, real_t* moment0, real_t* moment1);
    void ComputeBMomentGridBlock(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, m_ptr<const Field> fid_x);
};

class BMean : public BlockOperator {
   protected:
    lda_t  ida_;
    real_t sum_;

   public:
    explicit BMean() noexcept;
    explicit BMean(m_ptr<const Wavelet> interp) noexcept;

    void operator()(m_ptr<const ForestGrid> grid, m_ptr<const Field> fid_x, real_t* sum);
    void ComputeBMeanGridBlock(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, m_ptr<const Field> fid_x);
};

class BDiscreteMoment : public BlockOperator {
   protected:
    lda_t  ida_;
    real_t moment0_;
    real_t moment1_[3];

   public:
    explicit BDiscreteMoment() noexcept;
    explicit BDiscreteMoment(m_ptr<const Wavelet> interp) noexcept;

    void operator()(m_ptr<const ForestGrid> grid, m_ptr<const Field> fid_x, real_t* moment0, real_t* moment1);
    void ComputeBDiscreteMomentGridBlock(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, m_ptr<const Field> fid_x);
};

#endif