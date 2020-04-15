#include "subblock.hpp"

SubBlock::SubBlock(const lid_t gs, const lid_t stride, const lid_t start[3], const lid_t end[3]){
    gs_ = gs;
    stride_ = stride;
    for(int id=0; id<3; id++){
        end_[id] = end[id];
        start_[id] = start[id];
    }
}


void SubBlock::Reset(const lid_t gs, const lid_t stride, const lid_t start[3], const lid_t end[3]){
    gs_ = gs;
    stride_ = stride;
    for(int id=0; id<3; id++){
        end_[id] = end[id];
        start_[id] = start[id];
    }
}