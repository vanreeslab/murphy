#ifndef SRC_BLOCK_HPP_
#define SRC_BLOCK_HPP_

#include "field.hpp"
#include "murphy.hpp"

class Block {
   protected:
    sid_t level_;
    real_t xyz_[3]; 
    real_t hgrid_[3];

    datamap_t data_map_;

   public:
    Block(const real_t length, const real_t xyz[3], const sid_t level);
    ~Block();

    const real_t* xyz() const { return xyz_; }
    const real_t* hgrid() const { return hgrid_; }

    real_p       data(Field* fid);
    const real_p data(const Field* fid) const;

    void AddField(const qid_t* qid, Field* fid, nullptr_t ctx);
    void DeleteField(const qid_t* qid, Field* fid, nullptr_t ctx);
};

#endif  // SRC_BLOCK_HPP_