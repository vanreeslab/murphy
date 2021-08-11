#include "tools/toolsmpi.hpp"

#include "core/macros.hpp"

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

void GetRma(const level_t dlvl, const bidx_t shift[3],
            const MemLayout& layout_src, const MemSpan& span_src, MPI_Aint disp_src, rank_t src_rank,
            const MemLayout& layout_trg, const MemSpan& span_trg, MemData& data_trg, MPI_Win win) {
    m_assert(dlvl <= 1, "we cannot handle a difference in level > 1");
    m_assert(dlvl >= 0, "we cannot handle a level coarse ");
    m_assert(disp_src >= 0, "the displacement is not positive: %ld", disp_src);
    //-------------------------------------------------------------------------
    //................................................
    // get the corresponding MPI_Datatype for the target
    MPI_Datatype dtype_trg;
    ToMPIDatatype(layout_trg, span_trg, 1, &dtype_trg);

    //................................................
    // get the corresponding MPI_Datatype for the source
    MPI_Datatype dtype_src;
    ToMPIDatatype(layout_src, span_src, (bidx_t)pow(2, dlvl), &dtype_src);

    //................................................
    real_t*  local_trg = data_trg.ptr(span_trg.start[0], span_trg.start[1], span_trg.start[2]);
    MPI_Aint disp      = disp_src + layout_src.offset(span_src.start[0], span_src.start[1], span_src.start[2]);

    //..........................................................................
    int size_trg;
    MPI_Type_size(dtype_trg, &size_trg);
#ifndef NDEBUG
    int size_src;
    MPI_Type_size(dtype_src, &size_src);
    m_assert(size_src == size_trg, "the two sizes must be the same: %d vs %d", size_src, size_trg);
#endif
    // only perform the call if we expect something
    if (size_trg > 0) [[likely]] {
        MPI_Get(local_trg, 1, dtype_trg, src_rank, disp, 1, dtype_src, win);
    }

    // free the types
    MPI_Type_free(&dtype_trg);
    MPI_Type_free(&dtype_src);
    //-------------------------------------------------------------------------
}

/**
 * @brief use the MPI RMA Put function to copy the data from the disp_src to data_trg
 * 
 * @param dlvl the level gap = source level - target level (must be 0 or 1)
 * @param shift the shift, i.e. the position of the (0,0,0) of the target in the source framework (and resolution!)
 * @param block_src the memory layout corresponding to the source layout, only the ghost size and the stride are used
 * @param disp_src the displacement wrt to the target's window base pointer (the units are given by the disp_unit provided at the creation of the window on the target rank)
 * @param block_trg the memory layout corresponding to the target layout
 * @param data_trg the data memory pointer, i.e. the memory position of (0,0,0)
 * @param src_rank the rank of the source memory
 * @param win the window to use for the RMA calls
 */
void PutRma(const level_t dlvl, const bidx_t shift[3],
            const MemLayout& layout_src, const MemSpan& span_src, ConstMemData& data_src,
            const MemLayout& layout_trg, const MemSpan& span_trg, MPI_Aint disp_trg, rank_t trg_rank, MPI_Win win) {
    m_assert(dlvl <= 1, "we cannot handle a difference in level > 1");
    m_assert(dlvl >= 0, "we cannot handle a level coarse ");
    m_assert(disp_trg >= 0, "the displacement is not positive: %ld", disp_trg);
    //-------------------------------------------------------------------------
    //..........................................................................
    // get the corresponding MPI_Datatype for the target
    MPI_Datatype dtype_trg;
    ToMPIDatatype(layout_trg, span_trg, 1, &dtype_trg);

    //..........................................................................
    // get the corresponding MPI_Datatype for the source
    MPI_Datatype dtype_src;
    ToMPIDatatype(layout_src, span_src, (bidx_t)pow(2, dlvl), &dtype_src);

    //..........................................................................
    const real_t* local_src = data_src.ptr(span_src.start[0], span_src.start[1], span_src.start[2]);
    MPI_Aint      disp      = disp_trg + layout_trg.offset(span_trg.start[0], span_trg.start[1], span_trg.start[2]);

    //..........................................................................
    int size_trg;
    MPI_Type_size(dtype_trg, &size_trg);
#ifndef NDEBUG
    int size_src;
    MPI_Type_size(dtype_src, &size_src);
    m_assert(size_src == size_trg, "the two sizes must be the same: %d vs %d", size_src, size_trg);
#endif
    // only perform the call if we expect something
    if (size_trg > 0) [[likely]] {
        MPI_Put(local_src, 1, dtype_src, trg_rank, disp, 1, dtype_trg, win);
    }

    // free the types
    MPI_Type_free(&dtype_trg);
    MPI_Type_free(&dtype_src);
    //-------------------------------------------------------------------------
}