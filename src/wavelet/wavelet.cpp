#include "wavelet.hpp"

#include "core/forloop.hpp"

/**
 * @brief copy the data from data_src to data_trg
 * 
 * @warning we downsample the data if the levels do not match
 * @warning this is a wrapper to the @ref Wavelet::DoMagic_() function.
 * 
 * @param dlvl the difference of level: level_src - level_trg, only a difference of 0 or 1 is admissible
 * @param shift the position of the trg (0,0,0) in the src framework (and resolution!)
 * @param block_src description of data_src memory layout 
 * @param data_src the src data_ptr
 * @param block_trg descripiton of data_trg memory layout
 * @param data_trg the trg data_ptr
 */
void Wavelet::Copy(const level_t dlvl, const bidx_t shift[3], const MemLayout*  block_src, const_data_ptr data_src, const MemLayout*  block_trg, data_ptr data_trg) const {
    //-------------------------------------------------------------------------
    m_assert(dlvl == 0 || dlvl == 1, "only a difference of 0 or 1 is accepted, see the 2:1 constrain");
    // if not constant field, the target becomes its own constant field and the multiplication factor is 0.0
    DoMagic_(dlvl, true, shift, block_src, data_src, block_trg, data_trg, 0.0, nullptr);
    //-------------------------------------------------------------------------
}

/**
 * @brief interpolates (refine, coarsen or copy) the data from data_src to data_trg
 * 
 * This is a wrapper to the @ref Wavelet::DoMagic_() function.
 * The interpolation operation depends on the level difference.
 * 
 * @param dlvl the difference of level: level_src - level_trg, i.e. > 0 means coarsening, = 0 means copy and < 0 means refinement
 * @param shift the position of the trg (0,0,0) in the src framework (and resolution!)
 * @param block_src description of data_src memory layout 
 * @param data_src the 0-position of the src memory, i.e. the memory location of (0,0,0) for the source
 * @param block_trg descripiton of data_trg
 * @param data_trg the 0-position of the trg memory, i.e. the memory location of (0,0,0) for the target
 */
void Wavelet::Interpolate(const level_t dlvl, const bidx_t shift[3], const MemLayout*  block_src, const_data_ptr data_src, const MemLayout*  block_trg, data_ptr data_trg) const {
    //-------------------------------------------------------------------------
    // if not constant field, the target becomes its own constant field and the multiplication factor is 0.0
    DoMagic_(dlvl, false, shift, block_src, data_src, block_trg, data_trg, 0.0, nullptr);
    //-------------------------------------------------------------------------
}

/**
 * @brief interpolates the data from data_src and sum with the data_cst to data_trg: data_trg = alpha * data_cst + interp(data_src)
 * 
 * This is a wrapper to the @ref Wavelet::DoMagic_() function.
 * The interp() operation depends on the level difference.
 * 
 * @param dlvl the difference of level: level_src - level_trg, i.e. > 0 means coarsening, = 0 means copy and < 0 means refinement
 * @param shift the position of the trg (0,0,0) in the src framework (and resolution!)
 * @param block_src description of data_src memory layout 
 * @param data_src the 0-position of the src memory, i.e. the memory location of (0,0,0) for the source
 * @param block_trg descripiton of data_trg
 * @param data_trg the 0-position of the trg memory, i.e. the memory location of (0,0,0) for the target
 * @param alpha the multiplication factor for the constant memory
 * @param data_cst the 0-position of the constant memory, which follows the same layout as the target: block_trg
 * @param normal integers indicating the normal of the ghost layer. if not ghost, might be nullptr
 */
void Wavelet::Interpolate(const level_t dlvl, const bidx_t shift[3], const MemLayout*  block_src, const_data_ptr data_src, const MemLayout*  block_trg, data_ptr data_trg, const real_t alpha, const_data_ptr data_cst) const {
    //-------------------------------------------------------------------------
    // if not constant field, the target becomes its own constant field and the multiplication factor is 0.0
    DoMagic_(dlvl, false, shift, block_src, data_src, block_trg, data_trg, 0.0, nullptr);
    //-------------------------------------------------------------------------
}

/**
 * @brief interpolates the data from data_src and sum with the data_cst to data_trg: data_trg = alpha * data_cst + interp(data_src)
 * 
 * The interp() operation depends on the level difference.
 * This is a wrapper to the @ref Wavelet::DoMagic_() function.
 * 
 * @param dlvl the difference of level: level_src - level_trg, i.e. > 0 means coarsening, = 0 means copy and < 0 means refinement
 * @param force_copy if true, a copy is done instead of an interpolation, downsampling the data if necessary (then dlvl>0 is needed)
 * @param shift the position of the trg (0,0,0) in the src framework (and resolution!)
 * @param block_src description of data_src memory layout 
 * @param data_src the 0-position of the src memory, i.e. the memory location of (0,0,0) for the source
 * @param block_trg descripiton of data_trg
 * @param data_trg the 0-position of the trg memory, i.e. the memory location of (0,0,0) for the target
 * @param alpha the multiplication factor for the constant memory
 * @param data_cst the 0-position of the constant memory, which follows the same layout as the target: block_trg
 * @param normal integers indicating the normal of the ghost layer. if not ghost, might be nullptr
 */
void Wavelet::DoMagic_(const level_t dlvl, const bool force_copy, const bidx_t shift[3], const MemLayout*  block_src, const_data_ptr data_src, const MemLayout*  block_trg, data_ptr data_trg, const real_t alpha, const_data_ptr data_cst) const {
    m_assert(dlvl <= 1, "we cannot handle a difference in level > 1");
    m_assert(dlvl >= -1, "we cannot handle a level too coarse ");
    m_assert(!((alpha != 0.0 && data_cst.IsEmpty())), "if alpha = %e, the data_cst cannot be empty", alpha);
    //-------------------------------------------------------------------------
    // create the interpolation context
//     interp_ctx_t ctx;

//     // get memory details
//     for (sid_t id = 0; id < 3; id++) {
//         // the src starting and ending is place form the target point of view
// #ifndef NDEBUG
//         // for the target point of view, the start of the source is is in start - shift, same for the end
//         ctx.srcstart[id] = block_src->start(id) - shift[id];
//         ctx.srcend[id]   = block_src->end(id) - shift[id];
// #endif
//         ctx.trgstart[id] = block_trg->start(id);
//         ctx.trgend[id]   = block_trg->end(id);
//     }
//     ctx.srcstr = block_src->stride();
//     ctx.trgstr = block_trg->stride();

//     // store the multiplication factor and the normal
//     ctx.alpha = alpha;

//     // get the correct aligned arrays
//     // note: since the adresses refer to (0,0,0), we have a ghostsize of 0
//     ctx.sdata = data_src.Read(shift[0], shift[1], shift[2], 0, block_src->stride());
//     // ctx.cdata = data_cst;
//     ctx.tdata = data_trg;
    const_data_ptr sdata(data_src.Read(shift[0], shift[1], shift[2], 0), block_src->stride(), block_src->gs());

    bidx_t trgstart[3] = {block_trg->start(0), block_trg->start(1), block_trg->start(2)};
    bidx_t trgend[3]   = {block_trg->end(0), block_trg->end(1), block_trg->end(2)};

#ifndef NDEBUG
    bidx_t srcstart[3] = {block_src->start(0) - shift[0], block_src->start(1) - shift[1], block_src->start(2) - shift[2]};
    bidx_t srcend[3]   = {block_src->end(0) - shift[0], block_src->end(1) - shift[1], block_src->end(2) - shift[2]};

    interp_ctx_t ctx(sdata, block_src->stride(), srcstart, srcend, data_trg, block_trg->stride(), trgstart, trgend, alpha);
#else
    interp_ctx_t ctx(sdata, block_src->stride(), data_trg, block_trg->stride(), trgstart, trgend, alpha);
#endif
    

    // call the correct function
    if (dlvl == 0 || force_copy) {
        m_assert(dlvl >= 0, "the copy cannot be called with a negative dlvl argument, dlvl = %d, force copy = %d", dlvl, force_copy);
        Copy_(dlvl, &ctx);
    } else if (dlvl == -1) {
        RefineZeroDetails_(&ctx);
    } else if (dlvl > 0) {
        // m_log("src stride = %d, src start = %d %d %d, src end = %d %d %d, trg stride = %d , trg start = %d %d %d, trg end = %d %d %d, aalpha = %f", block_src->stride(), srcstart[0],srcstart[1],srcstart[2],
        //   srcend[0],srcend[1],srcend[2], block_trg->stride(),  trgstart[0], trgstart[1], trgstart[2], trgend[0], trgend[1], trgend[2], alpha);
        Coarsen_(&ctx);
    }
    //-------------------------------------------------------------------------
}

/**
 * @brief copy the value of the source memory to the target memory
 * 
 * @param dlvl the difference of level between the source and the target
 * @param ctx the interpolation context
 */
void Wavelet::Copy_(const level_t dlvl, const interp_ctx_t*  ctx) const {
    m_assert(dlvl <= 1, "we cannot handle a difference in level > 1");
    m_assert(dlvl >= 0, "we cannot handle a level coarse ");
    //-------------------------------------------------------------------------
    const bidx_t scaling = pow(2, dlvl);
    // const real_t alpha   = ctx->alpha;

    const bidx_t start[3] = {ctx->trgstart[0], ctx->trgstart[1], ctx->trgstart[2]};
    const bidx_t end[3]   = {ctx->trgend[0], ctx->trgend[1], ctx->trgend[2]};

    const bidx_t srcstr   = ctx->sdata.stride();
    const bidx_t trgstr   = ctx->tdata.stride();

    real_t*       tdata = ctx->tdata.Write();
    const real_t* sdata = ctx->sdata.Read();

    auto op = [=, &tdata](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
        // check we do not have to take the constant into account (not coded yet, weird behavior while saying trg = 0.0 * trg)
        // m_assert((ctx->cdata == nullptr), "the constant data must be nullptr for the moment");
        // check the accesses
        m_assert(((scaling * i0) >= ctx->srcstart[0]) && ((scaling * i0) < ctx->srcend[0]), "the source domain is too small in dir 0: %d >= %d and %d<%d", i0, ctx->srcstart[0], i0, ctx->srcend[0]);
        m_assert(((scaling * i1) >= ctx->srcstart[1]) && ((scaling * i1) < ctx->srcend[1]), "the source domain is too small in dir 1: %d >= %d and %d<%d", i1, ctx->srcstart[1], i1, ctx->srcend[1]);
        m_assert(((scaling * i2) >= ctx->srcstart[2]) && ((scaling * i2) < ctx->srcend[2]), "the source domain is too small in dir 2: %d >= %d and %d<%d", i2, ctx->srcstart[2], i2, ctx->srcend[2]);
        m_assert(((i0) >= ctx->trgstart[0]) && ((i0) < ctx->trgend[0]), "the target domain is too small in dir 0: %d >= %d and %d<%d", i0, ctx->trgstart[0], i0, ctx->trgend[0]);
        m_assert(((i1) >= ctx->trgstart[1]) && ((i1) < ctx->trgend[1]), "the target domain is too small in dir 1: %d >= %d and %d<%d", i1, ctx->trgstart[1], i1, ctx->trgend[1]);
        m_assert(((i2) >= ctx->trgstart[2]) && ((i2) < ctx->trgend[2]), "the target domain is too small in dir 2: %d >= %d and %d<%d", i2, ctx->trgstart[2], i2, ctx->trgend[2]);

        // get the current parent's data
        // const data_ptr lcdata = ctx->cdata + m_sidx(ik0, ik1, ik2, 0, ctx->trgstr);
        const real_t* lsdata = sdata + m_idx(scaling * i0, scaling * i1, scaling * i2, 0, srcstr);
        real_t*       ltdata = tdata + m_idx(i0, i1, i2, 0, trgstr);

        // do the simple copy
        ltdata[0] = lsdata[0];

        // non-nan checks
        m_assert(lsdata[0] == lsdata[0], "cannot be nan");
        // m_assert(lcdata[0] == lcdata[0], "cannot be nan");
        m_assert(ltdata[0] == ltdata[0], "cannot be nan");
    };

    for_loop(&op, start, end);
    //-------------------------------------------------------------------------
}

/**
 * @brief use the MPI RMA Get function to copy the data from the disp_src to data_trg
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
void Wavelet::GetRma(const level_t dlvl, const bidx_t shift[3], const MemLayout*  block_src, MPI_Aint disp_src, const MemLayout*  block_trg, data_ptr data_trg, rank_t src_rank, MPI_Win win) const {
    m_assert(dlvl <= 1, "we cannot handle a difference in level > 1");
    m_assert(dlvl >= 0, "we cannot handle a level coarse ");
    m_assert(disp_src >= 0, "the displacement is not positive: %ld", disp_src);
    //-------------------------------------------------------------------------
    //................................................
    // get the corresponding MPI_Datatype for the target
    const bidx_t  trg_start[3] = {block_trg->start(0), block_trg->start(1), block_trg->start(2)};
    const bidx_t  trg_end[3]   = {block_trg->end(0), block_trg->end(1), block_trg->end(2)};
    MPI_Datatype dtype_trg;
    ToMPIDatatype(trg_start, trg_end, block_trg->stride(), 1, &dtype_trg);

    //................................................
    // get the corresponding MPI_Datatype for the source
    const bidx_t  scale        = (bidx_t)pow(2, dlvl);
    const bidx_t  src_start[3] = {shift[0] + block_trg->start(0) * scale, shift[1] + block_trg->start(1) * scale, shift[2] + block_trg->start(2) * scale};
    const bidx_t  src_end[3]   = {shift[0] + block_trg->end(0) * scale, shift[1] + block_trg->end(1) * scale, shift[2] + block_trg->end(2) * scale};
    MPI_Datatype dtype_src;
    ToMPIDatatype(src_start, src_end, block_src->stride(), scale, &dtype_src);

    //................................................
    m_assert(data_trg.stride() == block_trg->stride(), "Stride of the layout is not the same as the stride in the memory.. This is baaaad");
    real_t*  local_trg = data_trg.Write(trg_start[0], trg_start[1], trg_start[2], 0);
    MPI_Aint disp      = disp_src + m_zeroidx(0, block_src) + m_idx(src_start[0], src_start[1], src_start[2], 0, block_src->stride());
    // //#pragma omp critical

    int size_trg;
    MPI_Type_size(dtype_trg, &size_trg);
#ifndef NDEBUG
    int size_src;
    MPI_Type_size(dtype_src, &size_src);
    m_assert(size_src == size_trg, "the two sizes must be the same: %d vs %d", size_src, size_trg);
#endif
    // only perform the call if we expect something
    if (size_trg > 0) {
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
void Wavelet::PutRma(const level_t dlvl, const bidx_t shift[3], const MemLayout*  block_src, const_data_ptr data_src, const MemLayout*  block_trg, MPI_Aint disp_trg, rank_t trg_rank, MPI_Win win) const {
    m_assert(dlvl <= 1, "we cannot handle a difference in level > 1");
    m_assert(dlvl >= 0, "we cannot handle a level coarse ");
    m_assert(disp_trg >= 0, "the displacement is not positive: %ld", disp_trg);
    //-------------------------------------------------------------------------
    //................................................
    // get the corresponding MPI_Datatype for the target
    const bidx_t  trg_start[3] = {block_trg->start(0), block_trg->start(1), block_trg->start(2)};
    const bidx_t  trg_end[3]   = {block_trg->end(0), block_trg->end(1), block_trg->end(2)};
    MPI_Datatype dtype_trg;
    ToMPIDatatype(trg_start, trg_end, block_trg->stride(), 1, &dtype_trg);

    //................................................
    // get the corresponding MPI_Datatype for the source
    const bidx_t  scale        = (bidx_t)pow(2, dlvl);
    const bidx_t  src_start[3] = {shift[0] + block_trg->start(0) * scale, shift[1] + block_trg->start(1) * scale, shift[2] + block_trg->start(2) * scale};
    const bidx_t  src_end[3]   = {shift[0] + block_trg->end(0) * scale, shift[1] + block_trg->end(1) * scale, shift[2] + block_trg->end(2) * scale};
    MPI_Datatype dtype_src;
    ToMPIDatatype(src_start, src_end, block_src->stride(), scale, &dtype_src);

    //................................................
    m_assert(data_src.stride() == block_src->stride(), "Stride of the layout is not the same as the stride in the memory.. This is baaaad");
    const real_t* local_src = data_src.Read(src_start[0], src_start[1], src_start[2], 0);
    MPI_Aint      disp      = disp_trg + m_zeroidx(0, block_trg) + m_idx(trg_start[0], trg_start[1], trg_start[2], 0, block_trg->stride());
    //#pragma omp critical
    int size_trg;
    MPI_Type_size(dtype_trg, &size_trg);
#ifndef NDEBUG
    int size_src;
    MPI_Type_size(dtype_src, &size_src);
    m_assert(size_src == size_trg, "the two sizes must be the same: %d vs %d", size_src, size_trg);
#endif
    // only perform the call if we expect something
    if (size_trg > 0) {
        MPI_Put(local_src, 1, dtype_src, trg_rank, disp, 1, dtype_trg, win);
    }

    // free the types
    MPI_Type_free(&dtype_trg);
    MPI_Type_free(&dtype_src);
    //-------------------------------------------------------------------------
}

/**
 * @brief return the biggest detail coefficient (infinite norm) as a refinement/coarsening criterion
 * 
 * The max detail is computed on an extension of the memory layout ( by criterion_shift_front and criterion_shift_back).
 * The reason for that comes from the coarsening of a block and that the scaling values are the function values
 * only if the contributing details are zero.
 *
 * @param block the block to analyze
 * @param data the data
 * @param smooth 
 * @return real_t the infinite norm of the max detail coefficient in the extended regio 
 */
real_t Wavelet::Criterion(/* source */ const MemLayout* const  block_src, const const_data_ptr& data,
                          /* target */ const MemLayout* const  detail_block) const {
    //-------------------------------------------------------------------------
    // get the detail coefficients
    real_t details_max = 0.0;
    
    Details(block_src, data, detail_block, data_ptr(nullptr), 0.0, &details_max);

    return details_max;
    //-------------------------------------------------------------------------
}

// /**
//  * @brief compute the Forward wavelet transform on a given block (scaling + detail) and return the max detail
//  *
//  * The details are computed on an extended domain as we need their values to drive the coarsening/refinement
//  *
//  * @param block_src the block for the data to be transformed (including ghosting points)
//  * @param block_trg the block for data to be computed
//  * @param data the data that will be transformed (in place)
//  * @param data_tmp the memory used as a copy (size @ref m_blockmemsize(1) )
//  * @return real_t
//  */
// real_t  Wavelet::FWT(const MemLayout* const  block, const data_ptr& data, const mem_ptr& tmp) const {
//     //-------------------------------------------------------------------------
//     data_ptr data_tmp = tmp(0,block);

//     // copy the block
//     ToTempMemory(block, data, data_tmp);

//     // get the extended memory layout
//     lid_t start[3];
//     lid_t end[3];
//     for (lda_t id = 0; id < 3; id++) {
//         start[id] = block->start(id) - criterion_shift_front();
//         end[id]   = block->end(id) + criterion_shift_back();
//     }
//     const SubBlock detail_block(block->gs(), block->stride(), start, end);

//     // get the detail coefficients
//     real_t details_max = 0.0;
//     Details(&detail_block, data, data_tmp, std::numeric_limits<real_t>::max(), &details_max);

//     // get now the scaling coefficients in the normal block
//     interp_ctx_t ctx;
//     for (sid_t id = 0; id < 3; id++) {
// #ifndef NDEBUG
//         ctx.srcstart[id] = block->start(id) - nghost_front();
//         ctx.srcend[id]   = block->end(id) + nghost_back();
// #endif
//         ctx.trgstart[id] = block->start(id);
//         ctx.trgend[id]   = block->end(id);
//     }
//     ctx.srcstr = block->stride();
//     ctx.trgstr = block->stride();
//     ctx.sdata = data_tmp; // we use data_tmp to get the scaling
//     ctx.tdata = data;

//     Coarsen_(&ctx)

//     return details_max;
//     //-------------------------------------------------------------------------
// }

// real_t Wavelet::CriterionAndSmooth(const MemLayout* const  block, const data_ptr& data, const mem_ptr& detail, const real_t tol) const {
//     //-------------------------------------------------------------------------
// // get the extended memory layout
// lid_t start[3];
// lid_t end[3];
// for (lda_t id = 0; id < 3; id++) {
//     start[id] = block->start(id) - criterion_shift_front();
//     end[id]   = block->end(id) + criterion_shift_back();
// }
// const SubBlock detail_block(block->gs(), block->stride(), start, end);

// // reset the detail array
// // memset(detail(),0,m_blockmemsize(1)*sizeof(real_t));
// data_ptr detail_data = detail(0,block);

// // get the detail coefficients
// real_t details_max = 0.0;
// Details(&detail_block, data, detail_data, tol, &details_max);

// // smooth them
// // Smooth(&detail_block,detail_data,block,data);

// return details_max;
//     //-------------------------------------------------------------------------
// }

/**
 * @brief Compute the detail coefficients are return the max infinite norm
 * 
 * @param block_src describes the source data layout
 * @param data the source data
 * @param detail_block describes the detail data layout
 * @param detail the detail data, if non-null the details are stored, cfr @ref Detail_ for storing policy
 * @param tol if > 0, used to decide wether a detail must be stored cfr @ref Detail_ for storing policy
 * @param details_max the max of the infinite norm on the details
 */
void Wavelet::Details(/* source */ const MemLayout* const  block_src, const const_data_ptr& data,
                      /* target */ const MemLayout* const  detail_block, const data_ptr& detail, const real_t tol,
                      /* output*/ real_t*  details_max) const {
    //-------------------------------------------------------------------------
    // get memory details
//     interp_ctx_t ctx;
//     for (lda_t id = 0; id < 3; id++) {
// #ifndef NDEBUG
//         ctx.srcstart[id] = block_src->start(id);
//         ctx.srcend[id]   = block_src->end(id);
// #endif
//         ctx.trgstart[id] = detail_block->start(id);
//         ctx.trgend[id]   = detail_block->end(id);
//     }
//     ctx.srcstr = block_src->stride();
//     ctx.sdata  = data;
//     ctx.trgstr = detail_block->stride();
//     ctx.tdata  = detail;
//     // we do not neeed to store
//     ctx.alpha = tol;


    bidx_t trgstart[3] = {detail_block->start(0), detail_block->start(1), detail_block->start(2)};
    bidx_t trgend[3]   = {detail_block->end(0), detail_block->end(1), detail_block->end(2)};

#ifndef NDEBUG
    bidx_t srcstart[3] = {block_src->start(0), block_src->start(1) , block_src->start(2)};
    bidx_t srcend[3]   = {block_src->end(0) , block_src->end(1) , block_src->end(2)};

    interp_ctx_t ctx(data, block_src->stride(), srcstart, srcend, detail, detail_block->stride(), trgstart, trgend, tol);
#else
    interp_ctx_t ctx(data, block_src->stride(), detail, detail_block->stride(), trgstart, trgend, tol);
#endif
    // we go for the two norm over the block
    Detail_(&ctx, details_max);
    //-------------------------------------------------------------------------
}

// void Wavelet::Scalings(const MemLayout* const  detail_block, const const_data_ptr& data, const data_ptr& detail, const real_t tol, real_t*  details_max) const {
//     //-------------------------------------------------------------------------
//     // get memory details
//     interp_ctx_t ctx;
//     for (lda_t id = 0; id < 3; id++) {
// #ifndef NDEBUG
//         ctx.srcstart[id] = 0 - nghost_front();
//         ctx.srcend[id]   = M_N + nghost_back();
// #endif
//         ctx.trgstart[id] = detail_block->start(id);
//         ctx.trgend[id]   = detail_block->end(id);
//     }
//     ctx.srcstr = detail_block->stride();
//     ctx.sdata  = data;
//     ctx.trgstr = detail_block->stride();
//     ctx.tdata  = detail;
//     // we do not neeed to store
//     ctx.alpha = tol;

//     // m_log("tol = %e", tol);
//     // we go for the two norm over the block
//     Detail_(&ctx, details_max);
//     //-------------------------------------------------------------------------
// }

/**
 * @brief Remove the details from the information if the mask is != 0.0
 * 
 * We use first compute the detail coefficients using the wavelet transform, then
 * we compute the inverse wavelet transform on the details to remove and correct the field
 * 
 * @param block_src the region used a source term for the detail computation
 * @param block_trg the region which need to be smoothed
 * @param data the data on which we operate
 * @param detail_block the region of details that need to be computed
 * @param detail_mask the mask on the details to compute (1.0 means we remove the detail)
 */
void Wavelet::SmoothOnMask(/* source */ const MemLayout* const  block_src,
                           /* target */ const MemLayout* const  block_trg, const data_ptr& data,
                           /* detail */ const MemLayout* const  detail_block, const data_ptr& detail_mask) const {
    //-------------------------------------------------------------------------
    //.........................................................................
    // get the details
    real_t trash = 0.0;
    Details(block_src, data, detail_block, detail_mask, -1.0, &trash);

    //.........................................................................
    // smooth them
    // get memory details
//     interp_ctx_t ctx;
//     for (lda_t id = 0; id < 3; id++) {
// #ifndef NDEBUG
//         // detail block is too small as we need more points that are 0.0 so use the full block_src... bad, I knooow
//         ctx.srcstart[id] = block_src->start(id);
//         ctx.srcend[id]   = block_src->end(id);
// #endif
//         ctx.trgstart[id] = block_trg->start(id);
//         ctx.trgend[id]   = block_trg->end(id);
//     }
//     ctx.srcstr = detail_block->stride();
//     ctx.sdata  = detail_mask;  // this are the details
//     ctx.trgstr = block_trg->stride();
//     ctx.tdata  = data;  // this is the block

    bidx_t trgstart[3] = {block_trg->start(0), block_trg->start(1), block_trg->start(2)};
    bidx_t trgend[3]   = {block_trg->end(0), block_trg->end(1), block_trg->end(2)};

#ifndef NDEBUG
    bidx_t srcstart[3] = {block_src->start(0), block_src->start(1), block_src->start(2)};
    bidx_t srcend[3]   = {block_src->end(0), block_src->end(1), block_src->end(2)};

    interp_ctx_t ctx(detail_mask, detail_block->stride(), srcstart, srcend, data, block_trg->stride(), trgstart, trgend, 0);
#else
    interp_ctx_t ctx(detail_mask, detail_block->stride(), data, block_trg->stride(), trgstart, trgend, 0);
#endif 
    // smooth them
    Smooth_(&ctx);
    //-------------------------------------------------------------------------
}

/**
 * @brief overwrites the data located at details coefficient's place to enforce a 0 detail computation
 * 
 * @warning this operation is quite brutal and should be use carefully as it doesn't conserve any moment property
 * 
 * @param block_src the region used as a source
 * @param block_trg the region where the details must be overwritten
 * @param data the modified data
 */
void Wavelet::OverwriteDetails(/* source */ const MemLayout* const  block_src,
                               /* target */ const MemLayout* const  block_trg, const data_ptr& data) const {
    //-------------------------------------------------------------------------
    //.........................................................................
//     interp_ctx_t ctx;
//     for (lda_t id = 0; id < 3; id++) {
// #ifndef NDEBUG
//         ctx.srcstart[id] = block_src->start(id);
//         ctx.srcend[id]   = block_src->end(id);
// #endif
//         ctx.trgstart[id] = block_trg->start(id);
//         ctx.trgend[id]   = block_trg->end(id);
//     }
//     // source and target data are the same!!
//     ctx.srcstr = block_trg->stride();
//     ctx.sdata  = data;
//     ctx.trgstr = block_trg->stride();
//     ctx.tdata  = data;

    bidx_t trgstart[3] = {block_trg->start(0), block_trg->start(1), block_trg->start(2)};
    bidx_t trgend[3]   = {block_trg->end(0), block_trg->end(1), block_trg->end(2)};
#ifndef NDEBUG
    bidx_t srcstart[3] = {block_src->start(0), block_src->start(1), block_src->start(2)};
    bidx_t srcend[3]   = {block_src->end(0), block_src->end(1), block_src->end(2)};

    interp_ctx_t ctx(data, block_trg->stride(), srcstart, srcend, data, block_trg->stride(), trgstart, trgend, 0);
#else
    interp_ctx_t ctx(data, block_trg->stride(), data, block_trg->stride(), trgstart, trgend, 0);
#endif
    // let's go
    OverwriteDetailsDualLifting_(&ctx);
    //-------------------------------------------------------------------------
}

/**
 * @brief Compute and store the details in the data_trg field
 */
void Wavelet::StoreDetails(/* source */ const MemLayout* const  block_src, const const_data_ptr& data,
                           /* target */ const MemLayout* const  block_detail, const data_ptr& detail) const {
    //-------------------------------------------------------------------------
    // get memory details
//     interp_ctx_t ctx;
//     for (lda_t id = 0; id < 3; id++) {
// #ifndef NDEBUG
//         ctx.srcstart[id] = block_src->start(id);
//         ctx.srcend[id]   = block_src->end(id);
// #endif
//         ctx.trgstart[id] = block_detail->start(id);
//         ctx.trgend[id]   = block_detail->end(id);
//     }
//     ctx.srcstr = block_src->stride();
//     ctx.trgstr = block_detail->stride();
//     ctx.sdata  = data;
//     ctx.tdata  = detail;
//     // set alpha to a huuge value, will store everything beneath it
//     ctx.alpha = std::numeric_limits<real_t>::max();

    bidx_t trgstart[3] = {block_detail->start(0), block_detail->start(1), block_detail->start(2)};
    bidx_t trgend[3]   = {block_detail->end(0), block_detail->end(1), block_detail->end(2)};
#ifndef NDEBUG
    bidx_t srcstart[3] = {block_src->start(0), block_src->start(1), block_src->start(2)};
    bidx_t srcend[3]   = {block_src->end(0), block_src->end(1), block_src->end(2)};
    
    interp_ctx_t ctx(data, block_src->stride(), srcstart, srcend, detail, block_detail->stride(), trgstart, trgend, std::numeric_limits<real_t>::max());
#else
    interp_ctx_t ctx(data, block_src->stride(), detail, block_detail->stride(), trgstart, trgend, std::numeric_limits<real_t>::max());
#endif
    // compute
    real_t detail_max = 0.0;  // -> will be discarded
    Detail_(&ctx, &detail_max);
    //-------------------------------------------------------------------------
}
