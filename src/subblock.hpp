#ifndef SRC_SUBBLOCK_HPP_
#define SRC_SUBBLOCK_HPP_

#include "murphy.hpp"

class SubBlock {
   protected:
    lid_t gs_;
    lid_t stride_;
    lid_t start_[3];
    lid_t range_[3];

    real_p data_;

   public:
    lid_t gs() const { return gs_; }
    lid_t stride() const { return stride_; }
    lid_t start(const int id) const { return start_[id]; }
    lid_t range(const int id) const { return range_[id]; }

    real_p data() { return data_; }
};

class GhostBlock : public SubBlock {
   protected:
    lid_t  dlvl_;
    lid_t  shift_[3];
    real_p data_src_;

   public:
    lid_t dlvl() const { return dlvl_; }
    lid_t shift(const int id) const { return shift_[id]; }

    lid_t* shift() { return shift_; }
    real_p data_src() { return data_src_; }
};

#endif  // SRC_GHOSTBLOCK_HPP_