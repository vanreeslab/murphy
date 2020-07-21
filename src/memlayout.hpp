#ifndef SRC_MEMLAYOUT_HPP_
#define SRC_MEMLAYOUT_HPP_

#include "mpi.h"
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

    virtual ~MemLayout(){};

    /**
     * @brief converts the current MemLayout to it's corresponding MPI Datatype vector and commit it
     * 
     * -> we change the lowerbound to account for the start index!
     * 
     * @return MPI_Datatype* 
     */
    virtual MPI_Datatype ToMPIDatatype() const {
        m_begin;
        //-------------------------------------------------------------------------
        // the number of blocks = the number of elem in y * number of elem in z
        int count = (end(1) - start(1)) * (end(2) - start(2));
        // the block length = the number of elem in x
        int blocklength = end(0) - start(0);

        // get the non-shifted type
        MPI_Datatype dtype_tmp;
        MPI_Type_vector(count, blocklength, stride(), M_MPI_REAL, &dtype_tmp);
        MPI_Type_commit(&dtype_tmp);

        // account for the offset in memory by changing the lower bound and keep the extent
        MPI_Datatype dtype;
        MPI_Aint     lb, extent;
        MPI_Type_get_extent(dtype_tmp, &lb, &extent);
        // update the lower bound
        lb = MPI_Aint_add(lb, start(0) * start(1) * start(2) * sizeof(real_t));
        // change the type and create a new one
        MPI_Type_create_resized(dtype_tmp, lb, extent, &dtype);
        MPI_Type_commit(&dtype);

        // delete the temp datatype
        MPI_Type_free(&dtype_tmp);

        // return the result
        return dtype;
        //-------------------------------------------------------------------------
        m_end;
    };
};

#endif  // SRC_MEMLAYOUT_HPP_
