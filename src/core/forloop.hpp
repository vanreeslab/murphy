#ifndef SRC_FORLOOP_HPP_
#define SRC_FORLOOP_HPP_

#include "core/macros.hpp"
#include "core/types.hpp"

template <bidx_t S0, bidx_t E0, bidx_t S1 = S0, bidx_t E1 = E0, bidx_t S2 = S0, bidx_t E2 = E0, class T>
inline void for_loop(T* core) {
    //-------------------------------------------------------------------------
    for (bidx_t i2 = S2; i2 < E2; ++i2) {
        for (bidx_t i1 = S1; i1 < E1; ++i1) {
            for (bidx_t i0 = S0; i0 < E0; ++i0) {
                (*core)(i0, i1, i2);
            }
        }
    }
    //-------------------------------------------------------------------------
};

template <class T>
inline void for_loop(T* core, const bidx_t s0, const bidx_t e0, const bidx_t s1, const bidx_t e1, const bidx_t s2, const bidx_t e2) {
    //-------------------------------------------------------------------------
    for (bidx_t i2 = s2; i2 < e2; ++i2) {
        for (bidx_t i1 = s1; i1 < e1; ++i1) {
            for (bidx_t i0 = s0; i0 < e0; ++i0) {
                (*core)(i0, i1, i2);
            }
        }
    }
    //-------------------------------------------------------------------------
};

template <class T>
inline void for_loop(T* core, const bidx_t start[3], const bidx_t end[3]) {
    //-------------------------------------------------------------------------
    for_loop(core, start[0], start[1], start[2], end[0], end[1], end[2]);
    //-------------------------------------------------------------------------
};

template <class T>
inline void for_loop(T* core, const bidx_t start, const bidx_t end) {
    //-------------------------------------------------------------------------
    for_loop(core, start, start, start, end, end, end);
    //-------------------------------------------------------------------------
};

#endif  // SRC_FORLOOP_HPP_