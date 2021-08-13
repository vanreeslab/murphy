#ifndef SRC_OPERATOR_XBLAS_HPP
#define SRC_OPERATOR_XBLAS_HPP

#include "core/macros.hpp"
#include "core/types.hpp"
#include "operator/blockoperator.hpp"
#include "grid/field.hpp"
#include "grid/forestgrid.hpp"
#include "grid/cartblock.hpp"

class BMax : public BlockOperator {
   public:
    explicit BMax() noexcept;
    explicit BMax(const bidx_t* ghost_len) noexcept;

    real_t operator()(const ForestGrid* grid, const Field& fid_x) const;
    void   ComputeBMaxGridBlock(/* in */ const qid_t* qid, const CartBlock* block, const Field* fid_x,
                              /* out */ real_t* max) const;
};

class BMinMax : public BlockOperator {
   public:
    explicit BMinMax() noexcept;
    explicit BMinMax(const bidx_t* ghost_len) noexcept;

    void operator()(const ForestGrid* grid, const Field& fid_x, real_t* min, real_t* max) const;
    void ComputeBMinMaxGridBlock(const qid_t* qid, const CartBlock* block, const Field& fid_x, const lda_t ida, real_t res[2]) const;
};

class BMoment : public BlockOperator {
   public:
    explicit BMoment() noexcept;
    explicit BMoment(const bidx_t* ghost_len) noexcept;

    void operator()(const ForestGrid* grid, const Field& fid_x, real_t* moment0, real_t* moment1) const;
    void ComputeBMomentGridBlock(const qid_t* qid, const CartBlock* block, const Field& fid_x, const lda_t ida, real_t moments[4]) const;
};

class BAvg : public BlockOperator {
   public:
    explicit BAvg() noexcept;
    explicit BAvg(const bidx_t* ghost_len) noexcept;

    void operator()(const ForestGrid* grid, const Field& fid_x, real_t* sum) const;
    void ComputeBAvgGridBlock(const qid_t* qid, const CartBlock* block, const Field& fid_x, const lda_t ida, real_t* sum) const;
};

// class BDiscreteMoment : public BlockOperator {
//    protected:
//     lda_t  ida_;
//     real_t moment0_;
//     real_t moment1_[3];

//    public:
//     explicit BDiscreteMoment() noexcept;
//     explicit BDiscreteMoment(const Wavelet*  interp) noexcept;

//     void operator()(const ForestGrid*  grid, const Field*  fid_x, real_t* moment0, real_t* moment1);
//     void ComputeBDiscreteMomentGridBlock(const qid_t*  qid, GridBlock*  block, const Field*  fid_x);
// };

#endif