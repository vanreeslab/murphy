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
    Block(const real_t length[3],const real_t xyz[3],const sid_t level);
    ~Block();

    void AddField(const qid_t* qid, Field* fid, nullptr_t ctx);
    void DeleteField(const qid_t* qid, Field* fid, nullptr_t ctx);
};

#endif  // SRC_BLOCK_HPP_