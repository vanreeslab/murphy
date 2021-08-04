#ifndef SRC_CORE_MEMLAYOUT_HPP_
#define SRC_CORE_MEMLAYOUT_HPP_

#include "mpi.h"
#include "core/macros.hpp"
#include "core/types.hpp"

#include <cstring>

/**
 * @brief descibe a fundamental 3D memory layout and a region of interest (start-end)
 * 
 * WARNING: we consider a cell-decentered information,
 * 
 * In 2D, a memory layout can be represented as:
 *  (with "x" meaning a ghost point and "o" a block point)
 * ```
 *                                 stride
 *          <--------------------------------------------------->
 *          +-------+----------------------------------+---------+
 *          | gs[0] |             core                 |  gs[1]  |
 *          |<----->|<-------------------------------->||<------->|
 *          +-------x----------------------------------+---------+
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
 *  |    (-gs[0],-gs[0])
 *  |
 *  +-------> x = memory access
 * ```
 * 
 */
class MemLayout {
   public:
    [[nodiscard]] virtual bidx_t stride() const               = 0;  //!< the stride in memory, i.e. = gs[0] + core + gs[1]
    [[nodiscard]] virtual bidx_t core() const                 = 0;  //!< the number of block points
    [[nodiscard]] virtual bidx_t gs() const                   = 0;  //!< the ghost point size (in front of the block )
    [[nodiscard]] virtual bidx_t start(const lda_t ida) const = 0;  //!< the starting point for the region of interest
    [[nodiscard]] virtual bidx_t end(const lda_t ida) const   = 0;  //!< the end point of the region of interest

    virtual ~MemLayout() = default;  //!< declare the constructor as virtual to ensure destruction
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
inline void ToMPIDatatype(const bidx_t start[3], const bidx_t end[3], const bidx_t stride, const bidx_t scale, MPI_Datatype* xyz_type) {
    m_begin;
    m_assert(scale == 1 || scale == 2, "the scale must be 1 or 2: here: %d", scale);
    m_assert(start[0] <= end[0], "the end = %d is smaller than the start = %d", end[0], start[0]);
    m_assert(start[1] <= end[1], "the end = %d is smaller than the start = %d", end[1], start[1]);
    m_assert(start[2] <= end[2], "the end = %d is smaller than the start = %d", end[2], start[2]);
    //-------------------------------------------------------------------------
    // get how much is one real (in bytes)
    MPI_Aint stride_x, trash_lb;
    MPI_Type_get_extent(M_MPI_REAL, &trash_lb, &stride_x);

    MPI_Datatype x_type, xy_type;
    //................................................
    // do x type as a simple vector
    bidx_t count_x = (end[0] - start[0]);
    m_assert(count_x >= 0, "we at least need to take 1 element");
    m_assert(count_x <= stride, "we cannot take more element than the stride");
    // MPI_Aint stride_x = sizeof(real_t);
    // MPI_Type_create_hvector(count_x / scale, 1, stride_x * scale, M_MPI_REAL, &x_type);
    MPI_Type_vector(count_x / scale, 1, scale, M_MPI_REAL, &x_type);
    //................................................
    // do y type
    bidx_t   count_y  = (end[1] - start[1]);
    MPI_Aint stride_y = stride_x * stride;
    m_assert(count_y >= 0, "we at least need to take 1 element");
    m_assert(count_y <= stride, "we cannot take more element than the stride");
    MPI_Type_create_hvector(count_y / scale, 1, (MPI_Aint)(stride_y * scale), x_type, &xy_type);

    //................................................
    // do z type
    bidx_t   count_z  = (end[2] - start[2]);
    MPI_Aint stride_z = stride_y * stride;
    m_assert(count_z >= 0, "we at least need to take 1 element");
    m_assert(count_z <= stride, "we cannot take more element than the stride");
    MPI_Type_create_hvector(count_z / scale, 1, (MPI_Aint)(stride_z * scale), xy_type, xyz_type);
    //................................................
    // finally commit the type so it's ready to use
    MPI_Type_commit(xyz_type);
    MPI_Type_free(&xy_type);
    MPI_Type_free(&x_type);
    //-------------------------------------------------------------------------
    m_end;
};

inline void ToMPIDatatype(const bidx_t start[3], const bidx_t end[3], const bidx_t stride, MPI_Datatype* xyz_type) {
    m_begin;
    // m_assert(scale == 1 || scale == 2, "the scale must be 1 or 2: here: %d", scale);
    // m_assert(start[0] < end[0], "the end = %d is smaller than the start = %d", end[0], start[0]);
    // m_assert(start[1] < end[1], "the end = %d is smaller than the start = %d", end[1], start[1]);
    // m_assert(start[2] < end[2], "the end = %d is smaller than the start = %d", end[2], start[2]);
    //--------------------------------------------------------------------------
    const bidx_t type_start[3]    = {0, 0, 0};  // the start index is 0,0,0 as we handle the start before the call!
    const bidx_t size[3]    = {stride, stride, stride};
    const bidx_t subsize[3] = {end[2] - start[2], end[1] - start[1], end[0] - start[0]};
    // m_log("subarray: %d %d %d @ %d %d %d from %d %d %d", subsize[0], subsize[1], subsize[2], start[0], start[1], start[2], size[0], size[1], size[2]);
    // commit the new type
    MPI_Type_create_subarray(3, size, subsize, type_start, MPI_ORDER_C, M_MPI_REAL, xyz_type);
    MPI_Type_commit(xyz_type);
    //--------------------------------------------------------------------------
    m_end;
};

#endif  // SRC_MEMLAYOUT_HPP_
