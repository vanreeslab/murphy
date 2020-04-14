#include "subblock.hpp"

SubBlock::SubBlock(const lid_t gs, const lid_t stride, const lid_t start[3], const lid_t range[3]){
    gs_ = gs;
    stride_ = stride;
    for(int id=0; id<3; id++){
        range_[id] = range[id];
        start_[id] = start[id];
    }
}