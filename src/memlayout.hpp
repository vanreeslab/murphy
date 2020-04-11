#ifndef SRC_MEMLAYOUT_HPP_
#define SRC_MEMLAYOUT_HPP_

#include "murphy.hpp"

/**
 * @brief descibes the most fundamental memory layout capacity
 * 
 */
class MemLayout {
   public:
    virtual lid_t gs() const                = 0;
    virtual lid_t stride() const            = 0;
    virtual lid_t start(const int id) const = 0;
    virtual lid_t range(const int id) const = 0;
};

#endif  // SRC_MEMLAYOUT_HPP_
