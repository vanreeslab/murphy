#ifndef SRC_CORE_MEMSPAN_HPP_
#define SRC_CORE_MEMSPAN_HPP_

#include "core/types.hpp"
#include "core/macros.hpp"

struct MemLayout {
    size_t gs;         //!< number of ghost points in front of the block
    size_t shift;      //!< the shift to apply to the fastest rotating dimension to have the 0 aligned
    size_t stride[2];  //!< the strides of the data: fastest rotating stride in [0], user stride in [1]

    size_t n_elem;  //!< number of elements in the memory layout (= minimal allocation size)

    explicit MemLayout() = delete;
    explicit MemLayout(const lda_t n_dim, const bidx_t n_gs_front, const bidx_t n_block, const bidx_t n_gs_back = -1) noexcept;
};

struct MemSpan {
    bidx_t start[3] = {0, 0, 0};  //!< starting index for the region of interest
    bidx_t end[3]   = {0, 0, 0};  //!< ending index for the region of interest

    explicit MemSpan() = delete;
    explicit MemSpan(const bidx_t in_start, const bidx_t in_end) noexcept;
    explicit MemSpan(const bidx_t in_start[3], const bidx_t in_end[3]) noexcept;
};

/**
 * @brief convert a pair (MemSpan , MemLayout) into a MPI datatype equivalent
 * 
 * @warning we are unable to take the initial displacement into account, so you need to point the memory to the start location while using it!
 * @note see exemple 4.13, page 123 of the MPI standard 3.1
 * 
 * @param span the MemSpan to convert
 * @param layout the MemLayout to convert
 * @param scale the scale coefficient (1 or 2). If 2, we take one points out of 2 in each direction
 * @param xyz_type the corresponding datattype
 */
void ToMPIDatatype(const MemSpan& span, const MemLayout& layout, const bidx_t scale, MPI_Datatype* xyz_type) {
    m_begin;
    m_assert(scale == 1 || scale == 2, "the scale must be 1 or 2: here: %d", scale);
    m_assert(span.start[0] <= span.end[0], "the end = %d is smaller than the start = %d", span.end[0], span.start[0]);
    m_assert(span.start[1] <= span.end[1], "the end = %d is smaller than the start = %d", span.end[1], span.start[1]);
    m_assert(span.start[2] <= span.end[2], "the end = %d is smaller than the start = %d", span.end[2], span.start[2]);
    //--------------------------------------------------------------------------
    // get how much is one real (in bytes)
    MPI_Aint stride_x = sizeof(real_t);
#ifndef NDEBUG
    {
        MPI_Aint stride_lb, trash_lb;
        MPI_Type_get_extent(M_MPI_REAL, &trash_lb, &stride_lb);
        m_assert(stride_x == stride_lb, "the two strides should be the same... I am confused here: %ld vs %ld", stride_lb, stride_x);
    }
#endif

    MPI_Datatype x_type, xy_type;
    //................................................
    // do x type as a simple vector
    bidx_t count_x = (span.end[0] - span.start[0]);
    m_assert(count_x >= 0, "we at least need to take 1 element");
    m_assert(count_x <= layout.stride[0], "we cannot take more element than the stride");
    MPI_Type_create_hvector(count_x / scale, 1, (MPI_Aint)(stride_x * scale), M_MPI_REAL, &x_type);
    //................................................
    // do y type
    bidx_t   count_y  = (span.end[1] - span.start[1]);
    MPI_Aint stride_y = stride_x * layout.stride[0];
    m_assert(count_y >= 0, "we at least need to take 1 element");
    m_assert(count_y <= layout.stride[1], "we cannot take more element than the stride");
    MPI_Type_create_hvector(count_y / scale, 1, (MPI_Aint)(stride_y * scale), x_type, &xy_type);

    //................................................
    // do z type
    bidx_t   count_z  = (span.end[2] - span.start[2]);
    MPI_Aint stride_z = stride_y * layout.stride[1];
    m_assert(count_z >= 0, "we at least need to take 1 element");
    m_assert(count_z <= layout.stride[1], "we cannot take more element than the stride");
    MPI_Type_create_hvector(count_z / scale, 1, (MPI_Aint)(stride_z * scale), xy_type, xyz_type);
    //................................................
    // finally commit the type so it's ready to use
    MPI_Type_commit(xyz_type);
    MPI_Type_free(&x_type);
    MPI_Type_free(&xy_type);

    //--------------------------------------------------------------------------
    m_end;
};

#endif