#include "patch.hpp"

#include "core/macros.hpp"
#include "core/types.hpp"

Patch::Patch(const real_t origin[3], const real_t length[3], const level_t level) {
    level_ = level;
#pragma unroll 3
    for (lda_t id = 0; id < 3; ++id) {
        origin_[id] = origin[id];
        length_[id] = length[id];
    }
}
