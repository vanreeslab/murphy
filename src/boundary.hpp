#ifndef SRC_BOUNDARY_HPP_
#define SRC_BOUNDARY_HPP_

#include "murphy.hpp"
#include "subblock.hpp"

/**
 * @brief possible positions of a Boundary condition
 * 
 */
typedef enum bcloc_t {
    M_XM, /**< X - */
    M_XP, /**< X + */
    M_YM, /**< Y - */
    M_YP, /**< Y + */
    M_ZM, /**< Z - */
    M_ZP  /**< Z + */
} bcloc_t;

typedef enum bctype_t{
    M_BC_DIR,
    M_BC_NEU,
    M_BC_ZERO,
    M_BC_EXTRAP
}bctype_t;

class Boundary {
   public:
    virtual void operator()(SubBlock* block) = 0;
};

#endif  // SRC_BOUNDARY_HPP_
