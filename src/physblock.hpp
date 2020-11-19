#ifndef SRC_PHYSBLOCK_HPP_
#define SRC_PHYSBLOCK_HPP_

#include "murphy.hpp"
#include "subblock.hpp"
#include "boundary.hpp"

/**
 * @brief Physical block: a @ref SubBlock that will be used to apply a physical boundary condition
 * 
 */
class PhysBlock : public SubBlock {
   protected:
    sid_t iface_;  //!< ID of the normal to the face
   public:
    sid_t iface() { return iface_; }
    // sid_t dir() { return iface_ / 2; }
    PhysBlock(const iface_t iface, const MemLayout* block, const lid_t nghost_front, const lid_t nghost_back);
};

#endif  // SRC_PHYSBLOCK_HPP_
