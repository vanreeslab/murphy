#ifndef SRC_PHYSBLOCK_HPP_
#define SRC_PHYSBLOCK_HPP_

#include "core/macros.hpp"
#include "core/types.hpp"
#include "ghostblock.hpp"

/**
 * @brief Physical block: a @ref SubBlock that will be used to apply a physical boundary condition
 * 
 */
class PhysBlock : public GhostBlock {
   protected:
    const iface_t iface_;  //!< ID of the normal to the face

   public:
    PhysBlock() = delete;
    PhysBlock(const iface_t iface, const MemLayout* block);
    // PhysBlock(const iface_t iface, const MemLayout* block, const lid_t nghost_front, const lid_t nghost_back);

    iface_t iface() { return iface_; }
};

#endif  // SRC_PHYSBLOCK_HPP_
