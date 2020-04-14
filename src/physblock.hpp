#ifndef SRC_PHYS_BLOCK_HPP_
#define SRC_PHYS_BLOCK_HPP_

#include "boundary.hpp"
#include "gridblock.hpp"
#include "memlayout.hpp"
#include "murphy.hpp"

/**
 * @brief Physical block: a @ref SubBlock that will be used to apply a physical boundary condition
 * 
 */
class PhysBlock : public SubBlock {
   protected:
    sid_t normal_sign_[3]; /**< the sign of the normal to the fphy*/

   public:
    sid_t normal_sign(const sid_t id) const { return normal_sign_[id]; }

    PhysBlock(const sid_t iface, GridBlock* block);
};

#endif  // SRC_PHYS_BLOCK_HPP_