#include "patch.hpp"

#include "murphy.hpp"

Patch::Patch(const real_t origin[3], const real_t length[3], const lid_t level) {
    level_ = level;
    for (int id = 0; id < 3; id++) {
        origin_[id] = origin[id];
        length_[id] = length[id];
    }
}
