#include "interpolator.hpp"

#include "subblock.hpp"

/**
 * @brief copy the data from data_src to data_trg
 * 
 * @warning we downsample the data if the levels do not match
 * this is a wrapper to the @ref InterpolatingWavelet::DoMagic_() function.
 * 
 * @param dlvl the difference of level: level_src - level_trg, only a difference of o or 1 is possible
 * @param shift the position of the trg (0,0,0) in the src framework (and resolution!)
 * @param block_src description of data_src memory layout 
 * @param data_src the 0-position of the src memory, i.e. the memory location of (0,0,0) for the source
 * @param block_trg descripiton of data_trg
 * @param data_trg the 0-position of the trg memory, i.e. the memory location of (0,0,0) for the target
 */
void InterpolatingWavelet::Copy(const level_t dlvl, const lid_t shift[3], const MemLayout* block_src, const data_ptr data_src, const MemLayout* block_trg, data_ptr data_trg) const{
    m_assert(dlvl == 0 || dlvl == 1, "only a difference of 0 or 1 is accepted, see the 2:1 constrain");
    // if not constant field, the target becomes its own constant field and the multiplication factor is 0.0
    DoMagic_(dlvl, true, shift, block_src, data_src, block_trg, data_trg, 0.0, data_trg);
}

/**
 * @brief interpolates (refine, coarsen or copy) the data from data_src to data_trg
 * 
 * This is a wrapper to the @ref InterpolatingWavelet::DoMagic_() function.
 * The interpolation operation depends on the level difference.
 * 
 * @param dlvl the difference of level: level_src - level_trg, i.e. > 0 means coarsening, = 0 means copy and < 0 means refinement
 * @param shift the position of the trg (0,0,0) in the src framework (and resolution!)
 * @param block_src description of data_src memory layout 
 * @param data_src the 0-position of the src memory, i.e. the memory location of (0,0,0) for the source
 * @param block_trg descripiton of data_trg
 * @param data_trg the 0-position of the trg memory, i.e. the memory location of (0,0,0) for the target
 */
void InterpolatingWavelet::Interpolate(const level_t dlvl, const lid_t shift[3], const MemLayout* block_src, const data_ptr data_src, const MemLayout* block_trg, data_ptr data_trg) const{
    // if not constant field, the target becomes its own constant field and the multiplication factor is 0.0
    DoMagic_(dlvl, false, shift, block_src, data_src, block_trg, data_trg, 0.0, data_trg);
}

/**
 * @brief interpolates the data from data_src and sum with the data_cst to data_trg: data_trg = alpha * data_cst + interp(data_src)
 * 
 * This is a wrapper to the @ref InterpolatingWavelet::DoMagic_() function.
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
void InterpolatingWavelet::Interpolate(const level_t dlvl, const lid_t shift[3], const MemLayout* block_src, const data_ptr data_src, const MemLayout* block_trg, data_ptr data_trg, const real_t alpha, const data_ptr data_cst) const{
    // if not constant field, the target becomes its own constant field and the multiplication factor is 0.0
    DoMagic_(dlvl, false, shift, block_src, data_src, block_trg, data_trg, 0.0, data_trg);
}

/**
 * @brief interpolates the data from data_src and sum with the data_cst to data_trg: data_trg = alpha * data_cst + interp(data_src)
 * 
 * The interp() operation depends on the level difference.
 * This is a wrapper to the @ref InterpolatingWavelet::DoMagic_() function.
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
void InterpolatingWavelet::DoMagic_(const level_t dlvl, const bool force_copy, const lid_t shift[3], const MemLayout* block_src, const data_ptr data_src, const MemLayout* block_trg, data_ptr data_trg, const real_t alpha, const data_ptr data_cst) const {
    m_assert(dlvl <= 1, "we cannot handle a difference in level > 1");
    m_assert(dlvl >= -1, "we cannot handle a level too coarse ");
    // m_assert(nghost_back() <= M_GS, "the number of BLOCK_GS is too low, should be at least %d, here %d", nghost_back(), M_GS);
    // m_assert(nghost_front() <= M_GS, "the number of BLOCK_GS is too low, should be at least %d, here %d", nghost_front(), M_GS);
    //-------------------------------------------------------------------------
    // create the interpolation context
    interp_ctx_t ctx;
    // m_log("-----------------");
    // m_log("entering interpolator with shift = %d %d %d", shift[0], shift[1], shift[2]);
    // m_log("entering interpolator with srcstart = %d %d %d", block_src->start(0), block_src->start(1), block_src->start(2));
    // m_log("entering interpolator with srcend = %d %d %d", block_src->end(0), block_src->end(1), block_src->end(2));
    // m_log("entering interpolator with trgstart = %d %d %d", block_trg->start(0), block_trg->start(1), block_trg->start(2));
    // m_log("entering interpolator with trgend = %d %d %d", block_trg->end(0), block_trg->end(1), block_trg->end(2));

    // get memory details
    for (sid_t id = 0; id < 3; id++) {
        // the src starting and ending is place form the target point of view
#ifndef NDEBUG
        // for the target point of view, the start of the source is is in start - shift, same for the end
        ctx.srcstart[id] = block_src->start(id) - shift[id];
        ctx.srcend[id]   = block_src->end(id) - shift[id];
#endif
        ctx.trgstart[id] = block_trg->start(id);
        ctx.trgend[id]   = block_trg->end(id);
    }
    ctx.srcstr = block_src->stride();
    ctx.trgstr = block_trg->stride();

    // store the multiplication factor and the normal
    ctx.alpha = alpha;
    // ctx.normal = normal;

    // get the correct aligned arrays
    // note: since the adresses refer to (0,0,0), we have a ghostsize of 0
    ctx.sdata = data_src + m_midx(shift[0], shift[1], shift[2], 0, block_src);
    ctx.cdata = data_cst;
    ctx.tdata = data_trg;

    // call the correct function
    if (dlvl == 0 || force_copy) {
        m_assert(dlvl >= 0, "the copy cannot be called with a negative dlvl argument, dlvl = %d, force copy = %d", dlvl, force_copy);
        Copy_(dlvl, &ctx);
    } else if (dlvl == -1) {
        Refine_(&ctx);
    } else if (dlvl > 0) {
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
void InterpolatingWavelet::Copy_(const level_t dlvl, const interp_ctx_t* ctx)const {
    m_assert(dlvl <= 1, "we cannot handle a difference in level > 1");
    m_assert(dlvl >= 0, "we cannot handle a level coarse ");
    //-------------------------------------------------------------------------
    // ensure alignment for target, constant and source
    // m_assume_aligned(ctx->tdata);
    // // m_assume_aligned(ctx->sdata);
    // m_assume_aligned(ctx->cdata);

    const lid_t  scaling = pow(2, dlvl);
    const real_t alpha   = ctx->alpha;

    // do the copy
    for (lid_t ik2 = ctx->trgstart[2]; ik2 < ctx->trgend[2]; ik2++) {
        for (lid_t ik1 = ctx->trgstart[1]; ik1 < ctx->trgend[1]; ik1++) {
            for (lid_t ik0 = ctx->trgstart[0]; ik0 < ctx->trgend[0]; ik0++) {
                m_assert(((scaling * ik0) >= ctx->srcstart[0]) && ((scaling * ik0) < ctx->srcend[0]), "the source domain is too small in dir 0: %d >= %d and %d<%d", ik0, ctx->srcstart[0], ik0, ctx->srcend[0]);
                m_assert(((scaling * ik1) >= ctx->srcstart[1]) && ((scaling * ik1) < ctx->srcend[1]), "the source domain is too small in dir 1: %d >= %d and %d<%d", ik1, ctx->srcstart[1], ik1, ctx->srcend[1]);
                m_assert(((scaling * ik2) >= ctx->srcstart[2]) && ((scaling * ik2) < ctx->srcend[2]), "the source domain is too small in dir 2: %d >= %d and %d<%d", ik2, ctx->srcstart[2], ik2, ctx->srcend[2]);
                // get the current parent's data
                data_ptr       ltdata = ctx->tdata + m_sidx(ik0, ik1, ik2, 0, ctx->trgstr);
                const data_ptr lcdata = ctx->cdata + m_sidx(ik0, ik1, ik2, 0, ctx->trgstr);
                const data_ptr lsdata = ctx->sdata + m_sidx(scaling * ik0, scaling * ik1, scaling * ik2, 0, ctx->srcstr);
                // do the simple copy
                m_assert(((ik0) >= ctx->trgstart[0]) && ((ik0) < ctx->trgend[0]), "the target domain is too small in dir 0: %d >= %d and %d<%d", ik0, ctx->trgstart[0], ik0, ctx->trgend[0]);
                m_assert(((ik1) >= ctx->trgstart[1]) && ((ik1) < ctx->trgend[1]), "the target domain is too small in dir 1: %d >= %d and %d<%d", ik1, ctx->trgstart[1], ik1, ctx->trgend[1]);
                m_assert(((ik2) >= ctx->trgstart[2]) && ((ik2) < ctx->trgend[2]), "the target domain is too small in dir 2: %d >= %d and %d<%d", ik2, ctx->trgstart[2], ik2, ctx->trgend[2]);
                ltdata[0] = alpha * lcdata[0] + lsdata[0];
            }
        }
    }
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
void InterpolatingWavelet::GetRma(const level_t dlvl, const lid_t shift[3], const MemLayout* block_src, MPI_Aint disp_src, const MemLayout* block_trg, data_ptr data_trg, rank_t src_rank, MPI_Win win) const {
    m_assert(dlvl <= 1, "we cannot handle a difference in level > 1");
    m_assert(dlvl >= 0, "we cannot handle a level coarse ");
    m_assert(disp_src >= 0, "the displacement is not positive: %ld", disp_src);
    //-------------------------------------------------------------------------
    // m_log("-----------------");
    // m_log("entering interpolator with shift = %d %d %d", shift[0], shift[1], shift[2]);
    // m_log("entering interpolator with srcstart = %d %d %d", block_src->start(0), block_src->start(1), block_src->start(2));
    // m_log("entering interpolator with srcend = %d %d %d", block_src->end(0), block_src->end(1), block_src->end(2));
    // m_log("entering interpolator with trgstart = %d %d %d", block_trg->start(0), block_trg->start(1), block_trg->start(2));
    // m_log("entering interpolator with trgend = %d %d %d", block_trg->end(0), block_trg->end(1), block_trg->end(2));

    //................................................
    // get the corresponding MPI_Datatype for the target
    const lid_t  trg_start[3] = {block_trg->start(0), block_trg->start(1), block_trg->start(2)};
    const lid_t  trg_end[3]   = {block_trg->end(0), block_trg->end(1), block_trg->end(2)};
    MPI_Datatype dtype_trg;
    ToMPIDatatype(trg_start, trg_end, block_trg->stride(), 1, &dtype_trg);

    //................................................
    // get the corresponding MPI_Datatype for the source
    const lid_t  scale        = (lid_t)pow(2, dlvl);
    const lid_t  src_start[3] = {shift[0] + block_trg->start(0) * scale, shift[1] + block_trg->start(1) * scale, shift[2] + block_trg->start(2) * scale};
    const lid_t  src_end[3]   = {shift[0] + block_trg->end(0) * scale, shift[1] + block_trg->end(1) * scale, shift[2] + block_trg->end(2) * scale};
    MPI_Datatype dtype_src;
    ToMPIDatatype(src_start, src_end, block_src->stride(), scale, &dtype_src);

    //................................................
#ifndef NDEBUG
    int size_trg, size_src;
    MPI_Type_size(dtype_src, &size_src);
    MPI_Type_size(dtype_trg, &size_trg);
    m_assert(size_trg == size_src, "the two sizes must match: src = %d vs trg = %d", size_src, size_trg);
#endif
    data_ptr local_trg = data_trg + m_midx(trg_start[0], trg_start[1], trg_start[2], 0, block_trg);
    MPI_Aint disp      = disp_src + m_zeroidx(0, block_src) + m_midx(src_start[0], src_start[1], src_start[2], 0, block_src);
    // //#pragma omp critical
    MPI_Get(local_trg, 1, dtype_trg, src_rank, disp, 1, dtype_src, win);

    // free the types
    MPI_Type_free(&dtype_trg);
    MPI_Type_free(&dtype_src);
    //-------------------------------------------------------------------------
}

void InterpolatingWavelet::PutRma(const level_t dlvl, const lid_t shift[3], const MemLayout* block_src, const data_ptr ptr_src, const MemLayout* block_trg, MPI_Aint disp_trg, rank_t trg_rank, MPI_Win win)const {
    m_assert(dlvl <= 1, "we cannot handle a difference in level > 1");
    m_assert(dlvl >= 0, "we cannot handle a level coarse ");
    m_assert(disp_trg >= 0, "the displacement is not positive: %ld", disp_trg);
    //-------------------------------------------------------------------------
    //................................................
    // get the corresponding MPI_Datatype for the target
    const lid_t  trg_start[3] = {block_trg->start(0), block_trg->start(1), block_trg->start(2)};
    const lid_t  trg_end[3]   = {block_trg->end(0), block_trg->end(1), block_trg->end(2)};
    MPI_Datatype dtype_trg;
    ToMPIDatatype(trg_start, trg_end, block_trg->stride(), 1, &dtype_trg);

    //................................................
    // get the corresponding MPI_Datatype for the source
    const lid_t  scale        = (lid_t)pow(2, dlvl);
    const lid_t  src_start[3] = {shift[0] + block_trg->start(0) * scale, shift[1] + block_trg->start(1) * scale, shift[2] + block_trg->start(2) * scale};
    const lid_t  src_end[3]   = {shift[0] + block_trg->end(0) * scale, shift[1] + block_trg->end(1) * scale, shift[2] + block_trg->end(2) * scale};
    MPI_Datatype dtype_src;
    ToMPIDatatype(src_start, src_end, block_src->stride(), scale, &dtype_src);

    //................................................
#ifndef NDEBUG
    int size_trg, size_src;
    MPI_Type_size(dtype_src, &size_src);
    MPI_Type_size(dtype_trg, &size_trg);
    // m_verb("src size = %d and the trg size = %d", size_src, size_trg);
    m_assert(size_trg == size_src, "the two sizes must match: src = %d vs trg = %d", size_src, size_trg);
#endif
    data_ptr local_src = ptr_src + m_zeroidx(0, block_src) + m_midx(src_start[0], src_start[1], src_start[2], 0, block_src);
    MPI_Aint disp      = disp_trg + m_zeroidx(0, block_trg) + m_midx(trg_start[0], trg_start[1], trg_start[2], 0, block_trg);
    //#pragma omp critical
    MPI_Put(local_src, 1, dtype_src, trg_rank, disp, 1, dtype_trg, win);

    // free the types
    MPI_Type_free(&dtype_trg);
    MPI_Type_free(&dtype_src);
    //-------------------------------------------------------------------------
}

/**
 * @brief return the biggest detail coefficient (infinite norm) as a refinement/coarsening criterion
 * 
 * The max/min is computed on an extension of the memory Layout (by lifting_len/2).
 * The reason for that comes from the coarsening of a block. 
 * While coarsening a block, we implicitly assume that all detail coefficients involved in the refinement are 0.0.
 * By doing so, we can coarsen the block, forget the detail coefficients and still retreive a perfect information.
 * Practically, it means that the lifting step, i.e. the contribution of the detail coefficients to my scaling coeff, is useless.
 * 
 * The lifting step is driven by the lifting coefficient, whose total length is (2*Nt-1).
 * In front of the block, I need to ensure that the details from my neighbor are 0.0, which means (2*Nt-1)/2 detail are zero
 * At the back of the block, I need to ensure (2*Nt-1)/2 -1 details are 0 (as the last point is a detail)
 *
 * @param block the block to analyze
 * @param data the data
 * @return real_t the infinite norm of the max detail coefficient in the extended region
 */
real_t InterpolatingWavelet::Criterion(MemLayout* block, data_ptr data,mem_ptr data_tmp) const {
    //-------------------------------------------------------------------------
    // get the extended memory layout
    const lid_t lift_len = (2 * Nt() - 1);
    lid_t       start[3];
    lid_t       end[3];
    for (lda_t id = 0; id < 3; id++) {
        start[id] = block->start(id) - m_max(lift_len / 2, 0);
        end[id]   = block->end(id) + m_max((lift_len / 2) - 1, 0);
    }
    SubBlock extended_block(block->gs(), block->stride(), start, end);

    // get the detail coefficients
    real_t details_max = 0.0;
    Details(&extended_block, data,data_tmp, &details_max);

    return details_max;
    //-------------------------------------------------------------------------
}

/**
 * @brief compute the max detail coefficients of a given MemLayout
 * 
 * @param block the block on which we computed
 * @param data the memory pointer to the point (0,0,0) of that block
 * @param data_tmp the temp memory of size CoarseStride()^3, see the ghosting
 * @param details_max an array of size 8 that will contain the detail coefficients: dx, dy, dz, dxy, dyz, dxz, dxyz, mean
 */
void InterpolatingWavelet::Details(MemLayout* block, data_ptr data_block, mem_ptr data_tmp, real_t* details_max) const {
    //-------------------------------------------------------------------------
    // get the subblock for the coarse
    const lid_t nghost_front = ncriterion_front() / 2;
    const lid_t nghost_back  = ncriterion_back() / 2;
    const lid_t stride       = nghost_front + M_HN + nghost_back;
    const lid_t start        = -nghost_front;
    const lid_t end          = M_HN + nghost_back;
    m_assert(stride <= CoarseStride(), "there is not enough space in the tmp memory");
    SubBlock coarse_block(nghost_front, stride, start, end);

    // reset the memory to 0.0
    memset(data_tmp, 0, CoarseMemSize());
    mem_ptr tmp = data_tmp + m_zeroidx(0, &coarse_block);

    // m_log("reading from %d to %d");

    // // downsample the information, ghost included
    // for (lda_t i2 = start; i2 < end; ++i2) {
    //     for (lda_t i1 = start; i1 < end; ++i1) {
    //         for (lda_t i0 = start; i0 < end; ++i0) {
    //             tmp[m_sidx(i0, i1, i2, 0, coarse_block.stride())] = data_block[m_sidx(i0 * 2, i1 * 2, i2 * 2, 0, block->stride())];
    //         }
    //     }
    // }
    // need to compute here the physical boundary conditions...

    // get memory details
    interp_ctx_t ctx;
    for (int id = 0; id < 3; id++) {
#ifndef NDEBUG
        ctx.srcstart[id] = -1;
        ctx.srcend[id]   = -2;
#endif
        ctx.trgstart[id] = block->start(id);
        ctx.trgend[id]   = block->end(id);
    }
    ctx.srcstr = coarse_block.stride();
    ctx.sdata  = tmp;
    ctx.trgstr = block->stride();
    ctx.tdata  = data_block;
    Detail_(&ctx, details_max);
    //-------------------------------------------------------------------------
}

/**
 * @brief coarsen the values of the source memory to gather them in the target memory.
 * 
 * @tparam N the number of vanishing moment
 * @tparam Nt the order of interpolation
 * @param ctx the interpolation context
 */
void InterpolatingWavelet::Coarsen_(const interp_ctx_t* ctx) const {
    //-------------------------------------------------------------------------
    // assure alignment for the target, the source, the constant and the temp data
    // m_assume_aligned(ctx->tdata);
    m_assume_aligned(ctx->sdata);
    // m_assume_aligned(ctx->cdata);

    const real_t  alpha  = ctx->alpha;
    const lid_t   ha_lim = ha_half_lim();
    const_mem_ptr ha     = ha_filter();

    for (lid_t ik2 = ctx->trgstart[2]; ik2 < ctx->trgend[2]; ++ik2) {
        for (lid_t ik1 = ctx->trgstart[1]; ik1 < ctx->trgend[1]; ++ik1) {
            for (lid_t ik0 = ctx->trgstart[0]; ik0 < ctx->trgend[0]; ++ik0) {
                // do some checks
                m_assert(((2 * ik0 - ha_lim) >= (ctx->srcstart[0])) && ((2 * ik0 + ha_lim) <= ctx->srcend[0]), "the source domain is too small in dir 0: %d >= %d and %d < %d", 2 * ik0 - ha_lim, ctx->srcstart[0], 2 * ik0 + ha_lim, ctx->srcend[0]);
                m_assert(((2 * ik1 - ha_lim) >= (ctx->srcstart[1])) && ((2 * ik1 + ha_lim) <= ctx->srcend[1]), "the source domain is too small in dir 1: %d >= %d and %d < %d", 2 * ik1 - ha_lim, ctx->srcstart[1], 2 * ik1 + ha_lim, ctx->srcend[1]);
                m_assert(((2 * ik2 - ha_lim) >= (ctx->srcstart[2])) && ((2 * ik2 + ha_lim) <= ctx->srcend[2]), "the source domain is too small in dir 2: %d >= %d and %d < %d", 2 * ik2 - ha_lim, ctx->srcstart[2], 2 * ik2 + ha_lim, ctx->srcend[2]);
                //get the local adress of the source, the target and the constant
                data_ptr       ltdata = ctx->tdata + m_sidx(ik0, ik1, ik2, 0, ctx->trgstr);
                const data_ptr lcdata = ctx->cdata + m_sidx(ik0, ik1, ik2, 0, ctx->trgstr);
                const data_ptr lsdata = ctx->sdata + m_sidx(2 * ik0, 2 * ik1, 2 * ik2, 0, ctx->srcstr);

                // add the constant
                ltdata[0] = alpha * lcdata[0];
                // apply the filter
                for (sid_t id2 = -ha_lim; id2 <= ha_lim; id2++) {
                    for (sid_t id1 = -ha_lim; id1 <= ha_lim; id1++) {
                        for (sid_t id0 = -ha_lim; id0 <= ha_lim; id0++) {
                            ltdata[0] += lsdata[m_sidx(id0, id1, id2, 0, ctx->srcstr)] * ha[id0] * ha[id1] * ha[id2];
                        }
                    }
                }
            }
        }
    }
    //-------------------------------------------------------------------------
}

/**
 * @brief refine the source memory to get the associated target memory information
 * 
 * Here, we assume that the detail coefficients are null.
 * Hence, the values of the function are the scaling coefficient and we simply apply the dual-lifting scheme to obtain the missing information
 * 
 * @tparam N the number of vanishing moment
 * @tparam Nt the order of interpolation
 * @param ctx the interpolation context
 */
void InterpolatingWavelet::Refine_(const interp_ctx_t* ctx) const {
    //-------------------------------------------------------------------------
    // assure alignment for the target, the source, the constant and the temp data
    const real_t  alpha  = ctx->alpha;
    const sid_t   gs_lim = gs_half_lim();
    const_mem_ptr gs     = gs_filter();

    const lid_t  start[3] = {ctx->trgstart[0], ctx->trgstart[1], ctx->trgstart[2]};
    const lid_t  end[3]   = {ctx->trgend[0], ctx->trgend[1], ctx->trgend[2]};
    const real_t one      = 1.0;

    // for each of the data for the needed target
    for (lid_t ik2 = start[2]; ik2 < end[2]; ++ik2) {
        for (lid_t ik1 = start[1]; ik1 < end[1]; ++ik1) {
            for (lid_t ik0 = start[0]; ik0 < end[0]; ++ik0) {
                // get 0 if odd, 1 if even (even if negative!!)
                const sid_t iy = m_sign(ik1) * (ik1 % 2);
                const sid_t ix = m_sign(ik0) * (ik0 % 2);
                const sid_t iz = m_sign(ik2) * (ik2 % 2);
                m_assert(ix == 0 || ix == 1, "this are the two possible values");
                m_assert(iy == 0 || iy == 1, "this are the two possible values");
                m_assert(iz == 0 || iz == 1, "this are the two possible values");

                // get the target location
                data_ptr       ltdata = ctx->tdata + m_sidx(ik0, ik1, ik2, 0, ctx->trgstr);
                const data_ptr lcdata = ctx->cdata + m_sidx(ik0, ik1, ik2, 0, ctx->trgstr);

                //get the local adress of the source, a bit more complicated to handle the negative numbers
                const lid_t    ik0_s  = (ik0 - ix) / 2;
                const lid_t    ik1_s  = (ik1 - iy) / 2;
                const lid_t    ik2_s  = (ik2 - iz) / 2;
                const data_ptr lsdata = ctx->sdata + m_sidx(ik0_s, ik1_s, ik2_s, 0, ctx->srcstr);
                m_assert((ik0_s * 2) <= ik0, "if not, we made something wrong...: source = %d, target = %d", ik0_s, ik0);
                m_assert((ik1_s * 2) <= ik1, "if not, we made something wrong...: source = %d, target = %d", ik1_s, ik1);
                m_assert((ik2_s * 2) <= ik2, "if not, we made something wrong...: source = %d, target = %d", ik2_s, ik2);

                // get the filter, depending on if I am odd or even
                const_mem_ptr gs_x         = (ix == 1) ? (gs) : (&one);
                const_mem_ptr gs_y         = (iy == 1) ? (gs) : (&one);
                const_mem_ptr gs_z         = (iz == 1) ? (gs) : (&one);
                const sid_t   lim_start[3] = {(gs_lim)*ix, (gs_lim)*iy, (gs_lim)*iz};
                const sid_t   lim_end[3]   = {(gs_lim + 1) * ix, (gs_lim + 1) * iy, (gs_lim + 1) * iz};

                m_assert(((ik0 / 2 - lim_start[0]) >= ctx->srcstart[0]) && ((ik0 / 2 + lim_end[0]) <= ctx->srcend[0]), "the source domain is too small in dir 0: %d >= %d and %d < %d", ik0 - gs_lim, ctx->srcstart[0], ik0 + gs_lim, ctx->srcend[0]);
                m_assert(((ik1 / 2 - lim_start[1]) >= ctx->srcstart[1]) && ((ik1 / 2 + lim_end[1]) <= ctx->srcend[1]), "the source domain is too small in dir 1: %d >= %d and %d < %d", ik1 - gs_lim, ctx->srcstart[1], ik1 + gs_lim, ctx->srcend[1]);
                m_assert(((ik2 / 2 - lim_start[2]) >= ctx->srcstart[2]) && ((ik2 / 2 + lim_end[2]) <= ctx->srcend[2]), "the source domain is too small in dir 2: %d >= %d and %d < %d", ik2 - gs_lim, ctx->srcstart[2], ik2 + gs_lim, ctx->srcend[2]);

                // add the constant array
                ltdata[m_sidx(0, 0, 0, 0, ctx->trgstr)] = alpha * lcdata[m_sidx(0, 0, 0, 0, ctx->trgstr)];

                // if one dim is even, id = 0, -> gs[0] = 1 and that's it
                // if one dim is odd, id = 1, -> we loop on gs, business as usual
                for (sid_t id2 = -lim_start[2]; id2 <= lim_end[2]; ++id2) {
                    for (sid_t id1 = -lim_start[1]; id1 <= lim_end[1]; ++id1) {
                        for (sid_t id0 = -lim_start[0]; id0 <= lim_end[0]; ++id0) {
                            const real_t fact = gs_x[id0] * gs_y[id1] * gs_z[id2];
                            ltdata[m_sidx(0, 0, 0, 0, ctx->trgstr)] += fact * lsdata[m_sidx(id0, id1, id2, 0, ctx->srcstr)];
                        }
                    }
                }
            }
        }
    }
    //-------------------------------------------------------------------------
}

/**
 * @brief gets the detail coefficients of the wavelet. This approximates the local slope of the data
 * 
 * @tparam order 
 * @param ctx only the trgdata information are used, the source is considered empty
 * @param details_inf_norm the maximum of the local detail coefficients
 */
void InterpolatingWavelet::Detail_(const interp_ctx_t* ctx, real_t* details_max) const {
    //-------------------------------------------------------------------------
    const sid_t   gs_lim   = gs_half_lim();
    const_mem_ptr gs       = gs_filter();
    const real_t  one      = 1.0;
    const lid_t   start[3] = {ctx->trgstart[0], ctx->trgstart[1], ctx->trgstart[2]};
    const lid_t   end[3]   = {ctx->trgend[0], ctx->trgend[1], ctx->trgend[2]};

    // for each of the data for the considered children
    (*details_max) = 0.0;
    for (lid_t ik2 = start[2]; ik2 < end[2]; ++ik2) {
        for (lid_t ik1 = start[1]; ik1 < end[1]; ++ik1) {
            for (lid_t ik0 = start[0]; ik0 < end[0]; ++ik0) {
            
                // get 0 if odd, 1 if even (even if negative!!)
                const lda_t iy = m_sign(ik1) * (ik1 % 2);
                const lda_t ix = m_sign(ik0) * (ik0 % 2);
                const lda_t iz = m_sign(ik2) * (ik2 % 2);
                m_assert(ix == 0 || ix == 1, "this are the two possible values");
                m_assert(iy == 0 || iy == 1, "this are the two possible values");
                m_assert(iz == 0 || iz == 1, "this are the two possible values");

                // get the nearest even data
                const lid_t ik0_s = (ik0 - ix);
                const lid_t ik1_s = (ik1 - iy);
                const lid_t ik2_s = (ik2 - iz);

                // get it's location
                const_mem_ptr ltdata = ctx->tdata + m_sidx(ik0_s, ik1_s, ik2_s, 0, ctx->trgstr);
                m_assume_aligned(ltdata);

                // get the filter, depending on if I am odd or even
                const_mem_ptr gs_x         = (ix == 1) ? (gs) : (&one);
                const_mem_ptr gs_y         = (iy == 1) ? (gs) : (&one);
                const_mem_ptr gs_z         = (iz == 1) ? (gs) : (&one);
                const sid_t   lim_start[3] = {(gs_lim)*ix, (gs_lim)*iy, (gs_lim)*iz};
                const sid_t   lim_end[3]   = {(gs_lim + 1) * ix, (gs_lim + 1) * iy, (gs_lim + 1) * iz};

                // if one dim is even, id = 0, -> gs[0] = 1 and that's it
                // if one dim is odd, id = 1, -> we loop on gs, business as usual
                real_t interp = 0.0;
                for (sid_t id2 = -lim_start[2]; id2 <= lim_end[2]; ++id2) {
                    for (sid_t id1 = -lim_start[1]; id1 <= lim_end[1]; ++id1) {
                        for (sid_t id0 = -lim_start[0]; id0 <= lim_end[0]; ++id0) {
                            const real_t fact = gs_x[id0] * gs_y[id1] * gs_z[id2];
                            interp += fact * ltdata[m_sidx(2 * id0, 2 * id1, 2 * id2, 0, ctx->trgstr)];
                        }
                    }
                }
                real_t detail = ctx->tdata[m_sidx(ik0, ik1, ik2, 0, ctx->trgstr)] - interp;

                // check the maximum
                (*details_max) = m_max(std::fabs(detail), (*details_max));
            }
        }
    }
    //-------------------------------------------------------------------------
}
