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

    // /**
    //  * @brief converts the current MemLayout to it's corresponding MPI Datatype vector and commit it
    //  * 
    //  * @warning we use here the POINTER view
    //  * 
    //  * -> we change the lowerbound to account for the start index!
    //  * 
    //  * @return MPI_Datatype* 
    //  */
    // virtual MPI_Datatype ptr_ToMPIDatatype() const {
    //     m_begin;
    //     //-------------------------------------------------------------------------
    //     // convert the 
    //     lid_t start_g[3] = {start(0) + gs(), start(1) + gs(), start(2) + gs()};
    //     lid_t end_g[3]   = {end(0) + gs(), end(1) + gs(), end(2) + gs()};

    //     MPI_Aint     lb, extent;
    //     MPI_Datatype x_type, xy_type, xyz_type;
    //     MPI_Datatype x_type_tmp, xy_type_tmp, xyz_type_tmp;

    //     //................................................
    //     // do x type
    //     int      count_x  = end_g[0] - start_g[0];
    //     MPI_Aint stride_x = sizeof(real_t);
    //     m_assert(count_x > 0, "we at least need to take 1 element");
    //     MPI_Type_create_hvector(count_x, 1, stride_x, M_MPI_REAL, &x_type_tmp);

    //     // shift to reach the correct spot
    //     MPI_Type_get_extent(x_type_tmp, &lb, &extent);
    //     lb = MPI_Aint_add(lb, stride_x * start_g[0]);
    //     MPI_Type_create_resized(x_type_tmp, lb, extent, &x_type);
    //     MPI_Type_commit(&x_type);

    //     //................................................
    //     // do y type
    //     int      count_y  = end_g[1] - start_g[1];
    //     MPI_Aint stride_y = sizeof(real_t) * stride();
    //     m_assert(count_y > 0, "we at least need to take 1 element");
    //     MPI_Type_create_hvector(count_y, 1, stride_y, x_type, &xy_type_tmp);

    //     // shift to reach the correct spot
    //     MPI_Type_get_extent(xy_type_tmp, &lb, &extent);
    //     lb = MPI_Aint_add(lb, stride_y * start_g[1]);
    //     MPI_Type_create_resized(xy_type_tmp, lb, extent, &xy_type);
    //     MPI_Type_commit(&xy_type);

    //     //................................................
    //     // do z type
    //     int      count_z  = end_g[2] - start_g[2];
    //     MPI_Aint stride_z = sizeof(real_t) * stride() * stride();
    //     m_assert(count_z > 0, "we at least need to take 1 element");
    //     MPI_Type_create_hvector(count_z, 1, stride_z, xy_type, &xyz_type_tmp);

    //     // shift to reach the correct spot
    //     MPI_Type_get_extent(xyz_type_tmp, &lb, &extent);
    //     lb = MPI_Aint_add(lb, stride_z * start_g[2]);
    //     MPI_Type_create_resized(xyz_type_tmp, lb, extent, &xyz_type);
    //     MPI_Type_commit(&xyz_type);

    //     MPI_Type_free(&x_type);
    //     MPI_Type_free(&xy_type);
    //     MPI_Type_free(&x_type_tmp);
    //     MPI_Type_free(&xy_type_tmp);
    //     MPI_Type_free(&xyz_type_tmp);


    //     return xyz_type_tmp;

    //     // //................................................
    //     // // do y type
    //     // int count_y = end_g[1] - start_g[1];
    //     // m_assert(count_y>0, "we at least need to take 1 element");
    //     // MPI_Type_create_hvector(count_y, 1, sizeof(real_t), x_type, &x_type);
    //     // MPI_Type_commit(&x_type);


    //     // // the number of blocks = the number of elem in y * number of elem in z
    //     // int count = (end_g[1] - start_g[1]) * (end_g[2] - start_g[2]);
    //     // // the block length = the number of elem in x
    //     // int blocklength = end_g[0] - start_g[0];

    //     // // get the non-shifted type
    //     // MPI_Datatype dtype_tmp;
    //     // MPI_Type_create_hvector(count, blocklength, stride() * sizeof(real_t), M_MPI_REAL, &dtype_tmp);
    //     // m_log("we take %d blocks of %d length with stride of %ld bytes ", count, blocklength, stride() * sizeof(real_t));
    //     // MPI_Type_commit(&dtype_tmp);

    //     // // account for the offset in memory by changing the lower bound and keep the extent
    //     // MPI_Datatype dtype;
    //     // MPI_Aint     lb, extent;

    //     // MPI_Type_get_extent(dtype_tmp, &lb, &extent);
    //     // // update the lower bound, works with global index as well
    //     // MPI_Aint shift_byte = m_sidx(start_g[0], start_g[1], start_g[2], 0, stride()) * sizeof(real_t);
    //     // m_log("lb = %ld and we add %ld (%d,%d,%d)", lb, shift_byte);
    //     // lb = MPI_Aint_add(lb, shift_byte);
    //     // // change the type and create a new one
    //     // MPI_Type_create_resized(dtype_tmp, lb, extent, &dtype);
    //     // MPI_Type_commit(&dtype);

    //     // // delete the temp datatype
    //     // MPI_Type_free(&dtype_tmp);

    //     // // return the result
    //     // return dtype;
    //     //-------------------------------------------------------------------------
    //     m_end;
    // };
};

#endif  // SRC_MEMLAYOUT_HPP_
