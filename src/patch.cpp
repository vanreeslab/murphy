#include "patch.hpp"

#include "murphy.hpp"

Patch::Patch(real_t origin[3], real_t length[3], lid_t level) {
    level_ = level;
    for (int id = 0; id < 3; id++) {
        origin_[id] = origin[id];
        length_[id] = length[id];
    }
}
