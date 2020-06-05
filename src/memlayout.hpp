#ifndef SRC_MEMLAYOUT_HPP_
#define SRC_MEMLAYOUT_HPP_

#include "murphy.hpp"

/**
 * @brief descibes the most fundamental 3D memory layout and a region of interest (start-end)
 * 
 * WARNING: we consider a cell-decentered information
 * 
 * In 2D, a memory layout can be represented as (with "x" meaning a ghost point and "o" a block point)
 * ```
 *                              stride
 *          <------------------------------------------------>
 *          +-------+----------------------------------+---------+ (stride-gs,stride-gs)
 *          |  gs   |                                  |    gs   |
 *          |<----->|                                  |<------->|
 *          +-------+----------------------------------+---------+
 *          |       |                                  |         |
 *          |       |                      end         |         |
 *          |       o    +------------------+          |         |
 *          |       |    |                  |          |         |
 *          |       o    |   region         |          |         |
 *          |       |    |      of          |          |         |
 *          |       o    |    interest      |          |         |
 *          |       |    |                  |          |         |
 *          |       o    +------------------+          |         |
 *          |       |  start                           |         |
 *          |       o    o    o    o    o    o    o    x    x    |
 *          |       |(0,0,0)                           |         |
 *          +-------o----o----o----o----o----o----o----x----x----|
 *          |       |                                  |         |
 *  y       |       |                                  |         |
 *  ^       +-------+----------------------------------+---------+
 *  |    (-gs,-gs)
 *  |
 *  +-------> x = memory access
 * ```
 * 
 */
class MemLayout {
   public:
    virtual lid_t gs() const                = 0;  //!< the ghost point size
    virtual lid_t stride() const            = 0;  //!< the stride in memory
    virtual lid_t start(const int id) const = 0;  //!< the starting point for the region of interest
    virtual lid_t end(const int id) const   = 0;  //!< the end point of the region of interest
};

#endif  // SRC_MEMLAYOUT_HPP_
