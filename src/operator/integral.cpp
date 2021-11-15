#include "integral.hpp"

#include "forloop.hpp"

template <>
void Integral<0>::ComputeIntegralBlock(const qid_t* qid, const CartBlock* block, const integrant_t op, real_t* result) const {
    //--------------------------------------------------------------------------
    // let's go!
    real_t lres = 0.0;
    auto   int0 = [=, &lres](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
        const real_t value = (op)(i0, i1, i2, block);
        lres += value;
    };
    // create a new span that will take the last GP
    for_loop(&int0, span_);
    const real_t vol = block->hgrid(0) * block->hgrid(1) * block->hgrid(2);
    (*result) += vol * lres;
    //--------------------------------------------------------------------------
}
template <>
void Integral<2>::ComputeIntegralBlock(const qid_t* qid, const CartBlock* block, const integrant_t op, real_t* result) const {
    //--------------------------------------------------------------------------
    // get the border information
    bidx_t dom_start[3] = {span_.start[0], span_.start[1], span_.start[2]};
    bidx_t dom_end[3]   = {span_.end[0], span_.end[1], span_.end[2]};

    // let's go!
    real_t lres = 0.0;
    auto   int2 = [=, &lres](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
        // if we are on the edge in 1 dim -> 1/2, 2 dims -> 1/4 and 3 dims -> 1/8
        const bool is_edge[3] = {(dom_start[0] == i0) || (i0 == dom_end[0]),
                                 (dom_start[1] == i1) || (i1 == dom_end[1]),
                                 (dom_start[2] == i2) || (i2 == dom_end[2])};

        const real_t fact  = 1.0 / ((1.0 + is_edge[0]) * (1.0 + is_edge[1]) * (1.0 + is_edge[2]));
        const real_t value = (op)(i0, i1, i2, block);
        lres += value * fact;
    };
    // create a new span that will take the last GP
    const bidx_t span_shift[2][3] = {{0, 0, 0}, {-1, -1, -1}};
    MemSpan      extended_span(&span_, span_shift);
    for_loop(&int2, extended_span);

    const real_t vol = block->hgrid(0) * block->hgrid(1) * block->hgrid(2);
    (*result) += lres * vol;
    //--------------------------------------------------------------------------
}
// template <>
// Integral<4>::ComputeIntegralBlock(const qid_t* qid, const CartBlock* block, const Field* field, const lda_t ida, const integrant_t* op) const;
