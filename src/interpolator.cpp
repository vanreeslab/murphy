#include "interpolator.hpp"

// void Interpolator::operator()(GhostBlock* block_trg) {
//     m_begin;
//     //-------------------------------------------------------------------------
//     Interpolator::operator()(block_trg->dlvl(), block_trg->shift(), block_src, block_trg);
//     //-------------------------------------------------------------------------
//     m_end;
// }

void Interpolator::interpolate_(const lid_t dlvl, const lid_t shift[3], MemLayout* block_src, real_p data_src, MemLayout* block_trg, real_p data_trg) {
    m_begin;
    m_assert(dlvl <= 2, "we cannot handle a difference in level > 2");
    m_assert(dlvl >= -1, "we cannot handle a level too coarse ");
    //-------------------------------------------------------------------------
    // get memory details
    for (int id = 0; id < 3; id++) {
        // the parent starting and ending is place form the child point of view
        srcstart_[id] = block_src->start(id) - shift[id];
        srcend_[id]   = block_src->start(id) + block_src->range(id) - shift[id];
        trgstart_[id] = block_trg->start(id);
        trgend_[id]   = block_trg->start(id) + block_trg->range(id);
    }
    srcgs_  = block_src->gs();
    trggs_  = block_trg->gs();
    srcstr_ = block_src->stride();
    trgstr_ = block_trg->stride();

    // get the correct aligned arrays
    real_p data_src_ = data_src + m_midx(shift[0], shift[1], shift[2], 0, block_src);
    real_p data_trg_ = data_trg + m_midx(0, 0, 0, 0, block_trg);

    // call the correct function
    if (dlvl == -1) {
        Refine_(data_src_, data_trg_);
    } else if (dlvl == 0) {
        Copy_(data_src_, data_trg_);
    } else if (dlvl > 0) {
        Coarsen_(dlvl, data_src_, data_trg_);
    }
    //-------------------------------------------------------------------------
    m_end;
}