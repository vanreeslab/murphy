#ifndef SRC_PHYSBLOCK_HPP_
#define SRC_PHYSBLOCK_HPP_

#include "gridblock.hpp"
#include "murphy.hpp"
#include "subblock.hpp"

/**
 * @brief Physical block: a @ref SubBlock that will be used to apply a physical boundary condition
 * 
 */
class PhysBlock : public SubBlock {
   protected:
    sid_t iface_;  //!< ID of the normal to the face
   public:
    sid_t iface() { return iface_; }
    sid_t dir() { return iface_ / 2; }
    PhysBlock(const sid_t iface, MemLayout* block, const sid_t nghost_front[3], const sid_t nghost_back[3]);
};

#endif  // SRC_PHYSBLOCK_HPP_
