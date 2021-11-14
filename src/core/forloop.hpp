#ifndef SRC_FORLOOP_HPP_
#define SRC_FORLOOP_HPP_

#include "core/macros.hpp"
#include "core/memspan.hpp"
#include "core/types.hpp"

template <typename T>
inline void for_loop(T*           kern,
                     const bidx_t s0, const bidx_t e0, const bidx_t j0,
                     const bidx_t s1, const bidx_t e1, const bidx_t j1,
                     const bidx_t s2, const bidx_t e2, const bidx_t j2) {
    //-------------------------------------------------------------------------
    const bidx_t n0   = ((e0 - s0) / j0);
    const bidx_t n01  = ((e1 - s1) / j1) * n0;
    const bidx_t n012 = ((e2 - s2) / j2) * n01;

    for (bidx_t i012 = 0; i012 < n012; ++i012) {
        const bidx_t i12 = (i012 % n01);
        const bidx_t i2  = s2 + (i012 / n01) * j2;
        const bidx_t i1  = s1 + (i12 / n0) * j1;
        const bidx_t i0  = s0 + (i12 % n0) * j0;

        (*kern)(i0, i1, i2);
    }
    //-------------------------------------------------------------------------
};

template <typename T>
inline void for_loop(T* kern, const MemSpan* const span) {
    //--------------------------------------------------------------------------
    for_loop(kern,
             span->start[0], span->end[0], 1,
             span->start[1], span->end[1], 1,
             span->start[2], span->end[2], 1);
    //--------------------------------------------------------------------------
};
template <typename T>
inline void for_loop(T* kern, const MemSpan* const span, const bidx_t jumps[3]) {
    //--------------------------------------------------------------------------
    for_loop(kern,
             span->start[0], span->end[0], jumps[0],
             span->start[1], span->end[1], jumps[1],
             span->start[2], span->end[2], jumps[2]);

    //--------------------------------------------------------------------------
};

template <typename T>
inline void for_loop(T* kern, const MemSpan& span) {
    //-------------------------------------------------------------------------
    for_loop(kern,
             span.start[0], span.end[0], 1,
             span.start[1], span.end[1], 1,
             span.start[2], span.end[2], 1);
    //-------------------------------------------------------------------------
};

template <typename T>
inline void for_loop(T* kern, const MemSpan& span, const bidx_t jumps[3]) {
    //-------------------------------------------------------------------------
    for_loop(kern,
             span.start[0], span.end[0], jumps[0],
             span.start[1], span.end[1], jumps[1],
             span.start[2], span.end[2], jumps[2]);
    //-------------------------------------------------------------------------
};

//==============================================================================
// this ones underneath should go!
//==============================================================================

template <bidx_t S0, bidx_t E0, bidx_t S1 = S0, bidx_t E1 = E0, bidx_t S2 = S0, bidx_t E2 = E0, class T>
inline void for_loop(T* core) {
    //-------------------------------------------------------------------------
    // for (bidx_t i2 = S2; i2 < E2; ++i2) {
    //     for (bidx_t i1 = S1; i1 < E1; ++i1) {
    //         for (bidx_t i0 = S0; i0 < E0; ++i0) {
    //             (*core)(i0, i1, i2);
    //         }
    //     }
    // }
    const bidx_t n0   = (E0 - S0);
    const bidx_t n01  = (E1 - S1) * n0;
    const bidx_t n012 = (E2 - S2) * n01;

    for (bidx_t i012 = 0; i012 < n012; ++i012) {
        const bidx_t i12 = (i012 % n01);
        const bidx_t i2  = S2 + (i012 / n01);
        const bidx_t i1  = S1 + (i12 / n0);
        const bidx_t i0  = S0 + (i12 % n0);

        (*core)(i0, i1, i2);
    }
    //-------------------------------------------------------------------------
};

// template <typn T>
// inline void for_loop(T* core, const bidx_t s0, const bidx_t e0, const bidx_t s1, const bidx_t e1, const bidx_t s2, const bidx_t e2) {
//     //-------------------------------------------------------------------------
//     // for (bidx_t i2 = s2; i2 < e2; ++i2) {
//     //     for (bidx_t i1 = s1; i1 < e1; ++i1) {
//     //         for (bidx_t i0 = s0; i0 < e0; ++i0) {
//     //             (*core)(i0, i1, i2);
//     //         }
//     //     }
//     // }
//     const bidx_t n0   = (e0 - s0);
//     const bidx_t n01  = (e1 - s1) * n0;
//     const bidx_t n012 = (e2 - s2) * n01;

//     for (bidx_t i012 = 0; i012 < n012; ++i012) {
//         const bidx_t i12 = (i012 % n01);
//         const bidx_t i2  = s2 + (i012 / n01);
//         const bidx_t i1  = s1 + (i12 / n0);
//         const bidx_t i0  = s0 + (i12 % n0);

//         (*core)(i0, i1, i2);
//     }
//     //-------------------------------------------------------------------------
// };

template <typename T>
inline void for_loop(T* core, const bidx_t start[3], const bidx_t end[3]) {
    //-------------------------------------------------------------------------
    for_loop(core, start[0], end[0], 1, start[1], end[1], 1, start[2], end[2], 1);
    //-------------------------------------------------------------------------
};

template <typename T>
inline void for_loop(T* core, const bidx_t start, const bidx_t end) {
    //-------------------------------------------------------------------------
    for_loop(core, start, end, 1, start, end, 1, start, end, 1);
    //-------------------------------------------------------------------------
};

#endif  // SRC_FORLOOP_HPP_