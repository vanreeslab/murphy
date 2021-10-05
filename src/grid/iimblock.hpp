#ifndef SRC_IIM_IIMBLOCK_HPP_
#define SRC_IIM_IIMBLOCK_HPP_

#include "grid/gridblock.hpp"

class IIMBlock: public GridBlock {
  public:

    IIMBlock(const real_t length, const real_t xyz[3], const sid_t level) : 
        GridBlock(length, xyz, level) {}
    
    void PrintIIMBlockOrigin() const { 
        printf("In IIMBLock, block origin is (%1.3f, %1.3f, %1.3f)\n", xyz_[0], xyz_[1], xyz_[2]);
    }
};

#endif