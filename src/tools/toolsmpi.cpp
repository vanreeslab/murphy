#include "toolsmpi.hpp"

#include "core/macros.hpp"

/**
 * @brief convert a pair (MemSpan , MemLayout) into a MPI datatype equivalent
 * 
 * @warning we are unable to take the initial displacement into account, so you need to point the memory to the start location while using it!
 * @note see exemple 4.13, page 123 of the MPI standard 3.1
 * 
 * @param span the MemSpan to convert
 * @param layout the MemLayout to convert
 * @param scale the scale coefficient (1 or 2). If 2, we take one points out of 2 in each direction, i.e. the stride is x2 and the number of elements is /2
 * @param xyz_type the corresponding datattype
 */
void ToMPIDatatype(const MemLayout* layout, const MemSpan* span, const bidx_t scale, MPI_Datatype* xyz_type) {
    m_begin;
    m_assert(scale == 1 || scale == 2, "the scale must be 1 or 2: here: %d", scale);
    m_assert(span->start[0] <= span->end[0], "the end = %d is smaller than the start = %d", span->end[0], span->start[0]);
    m_assert(span->start[1] <= span->end[1], "the end = %d is smaller than the start = %d", span->end[1], span->start[1]);
    m_assert(span->start[2] <= span->end[2], "the end = %d is smaller than the start = %d", span->end[2], span->start[2]);
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
    bidx_t count_x = (span->end[0] - span->start[0]);
    m_assert(count_x >= 0, "we at least need to take 1 element");
    m_assert(count_x <= layout->stride[0], "we cannot take more element than the stride");
    MPI_Type_create_hvector(count_x / scale, 1, (MPI_Aint)(stride_x * scale), M_MPI_REAL, &x_type);
    //................................................
    // do y type
    bidx_t   count_y  = (span->end[1] - span->start[1]);
    MPI_Aint stride_y = stride_x * layout->stride[0];
    m_assert(count_y >= 0, "we at least need to take 1 element");
    m_assert(count_y <= layout->stride[1], "we cannot take more element than the stride");
    MPI_Type_create_hvector(count_y / scale, 1, stride_y * scale, x_type, &xy_type);

    //................................................
    // do z type
    bidx_t   count_z  = (span->end[2] - span->start[2]);
    MPI_Aint stride_z = stride_y * layout->stride[1];
    m_assert(count_z >= 0, "we at least need to take 1 element");
    m_assert(count_z <= layout->stride[1], "we cannot take more element than the stride");
    MPI_Type_create_hvector(count_z / scale, 1, stride_z * scale, xy_type, xyz_type);
    //................................................
    // finally commit the type so it's ready to use
    MPI_Type_commit(xyz_type);
    MPI_Type_free(&x_type);
    MPI_Type_free(&xy_type);
    //--------------------------------------------------------------------------
    m_end;
}
/**
 * @brief Perform an RMA GET operation on a given MPI Window
 * 
 * The source layout is REMOTE and is copied to the LOCAL target layout
 * 
 * @param dlvl the difference of levels, can be 0 or 1. if 1, we downsample at the source
 * @param shift the shift, i.e. the position of the (0,0,0) of the target in the source framework (and resolution!)
 * @param layout_src the source MemLayout describing the source memory
 * @param span_src the source MemSpan
 * @param disp_src the displacement in the window where to find the source (it's equivalent to the raw pointer!)
 * @param src_rank the rank of the processor containing the source memory
 * @param layout_trg the target layout
 * @param span_trg the span of the target
 * @param data_trg the local data of the target
 * @param win the window that permits the RMA call
 */
void GetRma(const level_t dlvl, const bidx_t shift[3],
            const MemLayout* layout_src, const MemSpan* span_src, const MPI_Aint disp_src, rank_t src_rank,
            const MemLayout* layout_trg, const MemSpan* span_trg, const MemData* data_trg, MPI_Win win) {
    m_assert(dlvl <= 1, "we cannot handle a difference in level > 1");
    m_assert(dlvl >= 0, "we cannot handle a level coarse ");
    m_assert(disp_src >= 0, "the displacement is not positive: %ld", disp_src);
    //-------------------------------------------------------------------------
    //................................................
    // get the corresponding MPI_Datatype for the target
    MPI_Datatype dtype_trg;
    ToMPIDatatype(layout_trg, span_trg, (bidx_t)1, &dtype_trg);
    // get the local pointer for the local target
    real_t* local_trg = data_trg->ptr(span_trg->start[0], span_trg->start[1], span_trg->start[2]);

    //................................................
    // get the memory span corresponding to the data that is required (span_trg) on the source layout
    const bidx_t  scale        = (bidx_t)pow(2, dlvl);
    const lid_t   src_start[3] = {shift[0] + span_trg->start[0] * scale, shift[1] + span_trg->start[1] * scale, shift[2] + span_trg->start[2] * scale};
    const lid_t   src_end[3]   = {shift[0] + span_trg->end[0] * scale, shift[1] + span_trg->end[1] * scale, shift[2] + span_trg->end[2] * scale};
    const MemSpan shifted_span_src(src_start, src_end);

    m_assert(span_src->start[0] <= src_start[0] && src_start[0] <= span_src->end[0], "the src_start must be in the block: %d <= %d <= %d", span_src->start[0], src_start[0], span_src->end[0]);
    m_assert(span_src->start[1] <= src_start[1] && src_start[1] <= span_src->end[1], "the src_start must be in the block: %d <= %d <= %d", span_src->start[1], src_start[1], span_src->end[1]);
    m_assert(span_src->start[2] <= src_start[2] && src_start[2] <= span_src->end[2], "the src_start must be in the block: %d <= %d <= %d", span_src->start[2], src_start[2], span_src->end[2]);
    m_assert(span_src->start[0] <= src_end[0] && src_end[0] <= span_src->end[0], "the src_end must be in the block: %d <= %d <= %d", span_src->start[0], src_end[0], span_src->end[0]);
    m_assert(span_src->start[1] <= src_end[1] && src_end[1] <= span_src->end[1], "the src_end must be in the block: %d <= %d <= %d", span_src->start[1], src_end[1], span_src->end[1]);
    m_assert(span_src->start[2] <= src_end[2] && src_end[2] <= span_src->end[2], "the src_end must be in the block: %d <= %d <= %d", span_src->start[2], src_end[2], span_src->end[2]);

    // get the corresponding MPI_Datatype for the source (based on the shifted start/end!)
    // the shift is irrelevant for the type and the scale is handeled inside the function
    MPI_Datatype dtype_src;
    ToMPIDatatype(layout_src, &shifted_span_src, scale, &dtype_src);
    bidx_t   offset = layout_src->offset(src_start[0], src_start[1], src_start[2]);
    MPI_Aint disp   = disp_src + offset;
    m_assert(offset >= 0, "the offset = %d from %d %d %d must be >= 0", offset, src_start[0], src_start[1], src_start[2]);
    m_assert(offset < layout_src->n_elem, "the offset = %d from %d %d %d must be < %ld", offset, src_start[0], src_start[1], src_start[2], layout_src->n_elem);

    //..........................................................................
    int size_trg;
    MPI_Type_size(dtype_trg, &size_trg);
#ifndef NDEBUG
    int size_src;
    MPI_Type_size(dtype_src, &size_src);
    m_assert(size_src == size_trg, "the two sizes must be the same: %d vs %d, scale = %d", size_src, size_trg, scale);
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
 * @brief Perform an RMA PUT operation on a given MPI Window
 * 
 * The source layout is LOCAL and is copied to the REMOTE target layout
 * 
 * @param dlvl the difference of levels, can be 0 or 1. if 1, we downsample at the source
 * @param shift the shift, i.e. the position of the (0,0,0) of the target in the source framework (and resolution!)
 * @param layout_src 
 * @param span_src 
 * @param data_src 
 * @param layout_trg 
 * @param span_trg 
 * @param disp_trg 
 * @param trg_rank 
 * @param win 
 */
void PutRma(const level_t dlvl, const bidx_t shift[3],
            const MemLayout* layout_src, const MemSpan* span_src, const ConstMemData* data_src,
            const MemLayout* layout_trg, const MemSpan* span_trg, const MPI_Aint disp_trg, rank_t trg_rank, MPI_Win win) {
    m_assert(dlvl <= 1, "we cannot handle a difference in level > 1");
    m_assert(dlvl >= 0, "we cannot handle a level coarse ");
    m_assert(disp_trg >= 0, "the displacement is not positive: %ld", disp_trg);
    //-------------------------------------------------------------------------
    //..........................................................................
    // get the corresponding MPI_Datatype for the target
    MPI_Datatype dtype_trg;
    ToMPIDatatype(layout_trg, span_trg, 1, &dtype_trg);
    // get the displacement for the remote target
    MPI_Aint disp = disp_trg + layout_trg->offset(span_trg->start[0], span_trg->start[1], span_trg->start[2]);

    //..........................................................................
    const bidx_t scale = (bidx_t)pow(2, dlvl);
    // get the actual address for the local memory source
    const lid_t   src_start[3] = {shift[0] + span_trg->start[0] * scale, shift[1] + span_trg->start[1] * scale, shift[2] + span_trg->start[2] * scale};
    const lid_t   src_end[3]   = {shift[0] + span_trg->end[0] * scale, shift[1] + span_trg->end[1] * scale, shift[2] + span_trg->end[2] * scale};
    const MemSpan shifted_span_src(src_start, src_end);

    m_assert(span_src->start[0] <= src_start[0] && src_start[0] <= span_src->end[0], "the src_start must be in the block: %d <= %d <= %d", span_src->start[0], src_start[0], span_src->end[0]);
    m_assert(span_src->start[1] <= src_start[1] && src_start[1] <= span_src->end[1], "the src_start must be in the block: %d <= %d <= %d", span_src->start[1], src_start[1], span_src->end[1]);
    m_assert(span_src->start[2] <= src_start[2] && src_start[2] <= span_src->end[2], "the src_start must be in the block: %d <= %d <= %d", span_src->start[2], src_start[2], span_src->end[2]);
    m_assert(span_src->start[0] <= src_end[0] && src_end[0] <= span_src->end[0], "the src_end must be in the block: %d <= %d <= %d", span_src->start[0], src_end[0], span_src->end[0]);
    m_assert(span_src->start[1] <= src_end[1] && src_end[1] <= span_src->end[1], "the src_end must be in the block: %d <= %d <= %d", span_src->start[1], src_end[1], span_src->end[1]);
    m_assert(span_src->start[2] <= src_end[2] && src_end[2] <= span_src->end[2], "the src_end must be in the block: %d <= %d <= %d", span_src->start[2], src_end[2], span_src->end[2]);
    // get the corresponding MPI_Datatype for the source (based on the shifted span!)
    // the shift is irrelevant for the type and the scale is handeled inside the function
    MPI_Datatype dtype_src;
    ToMPIDatatype(layout_src, &shifted_span_src, scale, &dtype_src);
    // get the local source pointer
    const real_t* local_src = data_src->ptr(src_start[0], src_start[1], src_start[2]);

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