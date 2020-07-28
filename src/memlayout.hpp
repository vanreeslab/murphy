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
};

/**
 * @brief convert a 3D Memory layout into a MPI datatype equivalent
 * 
 * @warning we are unable to take the initial displacement into account, so you need to point the memory to the start location while using it
 * @note see exemple 4.13, page 123 of the MPI standard 3.1
 * 
 * @param start the starting point in memory (in a "data" layout, i.e. no ghost points included)
 * @param end the ending point in memory (in a "data" layout, i.e. no ghost points included)
 * @param gs the ghost size
 * @param stride the stide of the array
 * @param scale the scaling coefficient (1 or 2). If 2, we take one points out of 2 in each direction
 * @return MPI_Datatype the corresponding datattype
 */
inline void ToMPIDatatype(const lid_t start[3], const lid_t end[3], const lid_t gs, const lid_t stride, const lid_t scale, MPI_Datatype* xyz_type) {
    m_begin;
    m_assert(scale == 1 || scale == 2, "the scale must be 1 or 2: here: %d", scale);
    //-------------------------------------------------------------------------
    // convert the
    lid_t        start_g[3] = {(start[0] + gs), (start[1] + gs), (start[2] + gs)};
    lid_t        end_g[3]   = {(end[0] + gs), (end[1] + gs), (end[2] + gs)};
    MPI_Datatype x_type, xy_type;
    //................................................
    // do x type
    lid_t    count_x  = (end_g[0] - start_g[0]);
    MPI_Aint stride_x = sizeof(real_t);
    m_assert(count_x > 0, "we at least need to take 1 element");
    m_assert(count_x <= stride, "we cannot take more element than the stride");
    MPI_Type_create_hvector(count_x / scale, 1, stride_x * scale, M_MPI_REAL, &x_type);
    //................................................
    // do y type
    lid_t    count_y  = (end_g[1] - start_g[1]);
    MPI_Aint stride_y = stride_x * stride;
    m_assert(count_y > 0, "we at least need to take 1 element");
    m_assert(count_y <= stride, "we cannot take more element than the stride");
    MPI_Type_create_hvector(count_y / scale, 1, stride_y * scale, x_type, &xy_type);
    //................................................
    // do z type
    lid_t    count_z  = (end_g[2] - start_g[2]);
    MPI_Aint stride_z = stride_y * stride;
    m_assert(count_z > 0, "we at least need to take 1 element");
    m_assert(count_z <= stride, "we cannot take more element than the stride");
    MPI_Type_create_hvector(count_z / scale, 1, stride_z * scale, xy_type, xyz_type);
    MPI_Type_commit(xyz_type);
    //................................................
    // free the useless types and return
    MPI_Type_free(&x_type);
    MPI_Type_free(&xy_type);
    //-------------------------------------------------------------------------
    m_end;
};

#endif  // SRC_MEMLAYOUT_HPP_
