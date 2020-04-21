#ifndef SRC_PHYS_BLOCK_HPP_
#define SRC_PHYS_BLOCK_HPP_

#include "murphy.hpp"
#include "subblock.hpp"
#include "gridblock.hpp"

/**
 * @brief Physical block: a @ref SubBlock that will be used to apply a physical boundary condition
 * 
 */
class PhysBlock : public SubBlock {
   protected:
    sid_t iface_;  //!< ID of the normal to the face
    GridBlock* block_src_;  //!< store the reference to the underlying @ref GridBlock

   public:
    sid_t iface() { return iface_; };
    PhysBlock(const sid_t iface, GridBlock* block);
};

#endif  // SRC_PHYS_BLOCK_HPP_