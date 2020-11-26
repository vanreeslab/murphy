#ifndef SRC_FORLOOP_HPP_
#define SRC_FORLOOP_HPP_

#include "defs.hpp"

template <lid_t S0, lid_t E0, lid_t S1 = S0, lid_t E1 = E0, lid_t S2 = S0, lid_t E2 = E0, class T>
inline void for_loop(T* core) {
    //-------------------------------------------------------------------------
    for (lid_t i2 = S2; i2 < E2; ++i2) {
        for (lid_t i1 = S1; i1 < E1; ++i1) {
            for (lid_t i0 = S0; i0 < E0; ++i0) {
                (*core)(i0, i1, i2);
            }
        }
    }
    //-------------------------------------------------------------------------
};

#endif  // SRC_FORLOOP_HPP_