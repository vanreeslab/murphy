#ifndef SRC_PATCH_HPP_
#define SRC_PATCH_HPP_

#include "core/macros.hpp"
#include "core/types.hpp"

class Patch {
   protected:
    level_t level_;
    real_t  origin_[3] = {0.0, 0.0, 0.0};
    real_t  length_[3] = {0.0, 0.0, 0.0};

   public:
    explicit Patch(const real_t origin[3], const real_t length[3], const level_t level) noexcept;

    level_t level() const { return level_; }
    real_t* origin() { return origin_; }
    real_t* length() { return length_; }
    real_t  origin(const sid_t id) const { return origin_[id]; }
    real_t  length(const sid_t id) const { return length_[id]; }
};

#endif  // SRC_PATCH_HPP_
