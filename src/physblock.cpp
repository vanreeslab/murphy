#include "physblock.hpp"

PhysBlock::PhysBlock(const sid_t iface, GridBlock* block) {
    m_begin;
    //-------------------------------------------------------------------------
    gs_     = M_GS;
    stride_ = M_STRIDE;

    for (int id = 0; id < 3; id++) {
        normal_sign_[id] = 0;
        start_[id]       = -M_GS;
        range_[id]       = M_STRIDE;
    }

    // overwrite in the face direction
    const sid_t dir   = iface / 2;
    const sid_t sign  = iface % 2;  // sign = 1, -> we go plus, sign = 0 -> we go minus
    range_[dir]       = M_GS;
    start_[dir]       = (sign == 0) ? (-M_GS) : M_N;
    normal_sign_[dir] = (sign == 0) ? -1 : 1;
    //-------------------------------------------------------------------------
    m_end;
}