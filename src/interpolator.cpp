#include "interpolator.hpp"

/**
 * @brief
 * 
 * @note see exemple 4.13, page 123 of the MPI standard 3.1
 * 
 * @param start 
 * @param end 
 * @param gs 
 * @param stride 
 * @param scale 
 * @return MPI_Datatype 
 */
static inline MPI_Datatype ToMPIDatatype(const lid_t start[3], const lid_t end[3], const lid_t gs, const lid_t stride, const lid_t scale) {
    m_begin;
    m_assert(scale == 1 || scale == 2, "the scale must be 1 or 2: here: %d", scale);
    //-------------------------------------------------------------------------
    // convert the
    lid_t        start_g[3] = {(start[0] + gs), (start[1] + gs), (start[2] + gs)};
    lid_t        end_g[3]   = {(end[0] + gs), (end[1] + gs), (end[2] + gs)};
    MPI_Datatype x_type, xy_type, xyz_type;
    //................................................
    // do x type
    int      count_x  = (end_g[0] - start_g[0]);
    MPI_Aint stride_x = sizeof(real_t);
    m_assert(count_x > 0, "we at least need to take 1 element");
    m_assert(count_x <= stride, "we cannot take more element than the stride");
    m_log("we take %d elems with a stide of %ld", count_x / scale, stride_x * scale);
    MPI_Type_create_hvector(count_x / scale, 1, stride_x * scale, M_MPI_REAL, &x_type);
    //................................................
    // do y type
    int      count_y  = (end_g[1] - start_g[1]);
    MPI_Aint stride_y = stride_x * stride;
    m_assert(count_y > 0, "we at least need to take 1 element");
    m_assert(count_y <= stride, "we cannot take more element than the stride");
    m_log("we take %d elems with a stide of %ld", count_y / scale, stride_y * scale);
    MPI_Type_create_hvector(count_y / scale, 1, stride_y * scale, x_type, &xy_type);
    //................................................
    // do z type
    int      count_z  = (end_g[2] - start_g[2]);
    MPI_Aint stride_z = stride_y * stride;
    m_assert(count_z > 0, "we at least need to take 1 element");
    m_assert(count_z <= stride, "we cannot take more element than the stride");
    m_log("we take %d elems with a stide of %ld", count_z / scale, stride_z * scale);
    MPI_Type_create_hvector(count_z / scale, 1, stride_z * scale, xy_type, &xyz_type);
    MPI_Type_commit(&xyz_type);
    //................................................
    // free the useless types and return
    MPI_Type_free(&x_type);
    MPI_Type_free(&xy_type);
    return xyz_type;
    //-------------------------------------------------------------------------
    m_end;
};

/**
 * @brief copy the data from data_src to data_trg but downsampling the data if the levels do not match
 * 
 * @param dlvl the difference of level: level_src - level_trg, only a a positive difference is possible
 * @param shift the position of the trg (0,0,0) in the src framework (and resolution!)
 * @param block_src description of data_src memory layout 
 * @param data_src the 0-position of the src memory, i.e. the memory location of (0,0,0) for the source
 * @param block_trg descripiton of data_trg
 * @param data_trg the 0-position of the trg memory, i.e. the memory location of (0,0,0) for the target
 */
void Interpolator::Copy(const sid_t dlvl, const lid_t shift[3], MemLayout *block_src, real_p data_src, MemLayout *block_trg, real_p data_trg) {
    // if not constant field, the target becomes its own constant field and the multiplication factor is 0.0
    DoMagic_(dlvl, true, shift, block_src, data_src, block_trg, data_trg, 0.0, data_trg);
}

/**
 * @brief interpolates the data from data_src to data_trg
 * 
 * @param dlvl the difference of level: level_src - level_trg, i.e. > 0 means coarsening, = 0 means copy and < 0 means refinement
 * @param shift the position of the trg (0,0,0) in the src framework (and resolution!)
 * @param block_src description of data_src memory layout 
 * @param data_src the 0-position of the src memory, i.e. the memory location of (0,0,0) for the source
 * @param block_trg descripiton of data_trg
 * @param data_trg the 0-position of the trg memory, i.e. the memory location of (0,0,0) for the target
 */
void Interpolator::Interpolate(const sid_t dlvl, const lid_t shift[3], MemLayout *block_src, real_p data_src, MemLayout *block_trg, real_p data_trg) {
    // if not constant field, the target becomes its own constant field and the multiplication factor is 0.0
    DoMagic_(dlvl, false, shift, block_src, data_src, block_trg, data_trg, 0.0, data_trg);
}

/**
 * @brief interpolates the data from data_src and sum with the data_cst to data_trg: data_trg = alpha * data_cst + interp(data_src)
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
void Interpolator::Interpolate(const sid_t dlvl, const lid_t shift[3], MemLayout *block_src, real_p data_src, MemLayout *block_trg, real_p data_trg, const real_t alpha, real_p data_cst) {
    // if not constant field, the target becomes its own constant field and the multiplication factor is 0.0
    DoMagic_(dlvl, false, shift, block_src, data_src, block_trg, data_trg, 0.0, data_trg);
}

/**
 * @brief interpolates the data from data_src and sum with the data_cst to data_trg: data_trg = alpha * data_cst + interp(data_src)
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
void Interpolator::DoMagic_(const sid_t dlvl,const bool force_copy, const lid_t shift[3], MemLayout *block_src, real_p data_src, MemLayout *block_trg, real_p data_trg, const real_t alpha, real_p data_cst)
{
    m_assert(dlvl <= 1, "we cannot handle a difference in level > 1");
    m_assert(dlvl >= -1, "we cannot handle a level too coarse ");
    // m_assert(!(normal != nullptr && dlvl == 1), "if a normal is given, we only copy or refine");
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
    for (int id = 0; id < 3; id++) {
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
 * @param 
 * @param ctx the interpolation context
 */
void Interpolator::Copy_(const sid_t dlvl, const interp_ctx_t* ctx) {
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
                real_p       ltdata = ctx->tdata + m_sidx(ik0, ik1, ik2, 0, ctx->trgstr);
                const real_p lcdata = ctx->cdata + m_sidx(ik0, ik1, ik2, 0, ctx->trgstr);
                const real_p lsdata = ctx->sdata + m_sidx(scaling * ik0, scaling * ik1, scaling * ik2, 0, ctx->srcstr);
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
 * @brief copy the data from data_src to data_trg but downsampling the data if the levels do not match using MPI RMA calls
 * 
 * @param dlvl the difference of level: level_src - level_trg, only a a positive difference is possible
 * @param shift the position of the trg (0,0,0) in the src framework (and resolution!)
 * @param block_src description of data_src memory layout 
 * @param data_src the 0-position of the src memory, i.e. the memory location of (0,0,0) for the source
 * @param block_trg descripiton of data_trg
 * @param data_trg the 0-position of the trg memory, i.e. the memory location of (0,0,0) for the target
 */
void Interpolator::GetRma(const sid_t dlvl, const lid_t shift[3], MemLayout *block_src, MPI_Aint disp_src, MemLayout *block_trg, real_t *data_trg, int src_rank, MPI_Win win) {
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
    MPI_Datatype dtype_trg    = ToMPIDatatype(trg_start, trg_end, block_trg->gs(), block_trg->stride(), 1);
    // m_verb("the trg = %d %d %d to %d %d %d", trg_start[0], trg_start[1], trg_start[2], trg_end[0], trg_end[1], trg_end[2]);

    //................................................
    // get the corresponding MPI_Datatype for the target
    const lid_t  scale        = (lid_t)pow(2, dlvl);
    const lid_t  src_start[3] = {shift[0] + block_trg->start(0) * scale, shift[1] + block_trg->start(1) * scale, shift[2] + block_trg->start(2) * scale};
    const lid_t  src_end[3]   = {shift[0] + block_trg->end(0) * scale, shift[1] + block_trg->end(1) * scale, shift[2] + block_trg->end(2) * scale};
    MPI_Datatype dtype_src    = ToMPIDatatype(src_start, src_end, block_src->gs(), block_src->stride(), scale);
    // m_verb("the src = %d %d %d to %d %d %d", src_start[0], src_start[1], src_start[2], src_end[0], src_end[1], src_end[2]);

    //................................................
#ifndef NDEBUG
    int size_trg, size_src;
    MPI_Type_size(dtype_src, &size_src);
    MPI_Type_size(dtype_trg, &size_trg);
    // m_verb("src size = %d and the trg size = %d", size_src, size_trg);
    m_assert(size_trg == size_src, "the two sizes must match: src = %d vs trg = %d", size_src, size_trg);
#endif
    real_t * local_trg = data_trg + m_midx(trg_start[0], trg_start[1], trg_start[2], 0, block_trg);
    MPI_Aint disp      = disp_src + m_zeroidx(0, block_src) + m_midx(src_start[0], src_start[1], src_start[2], 0, block_src);
#pragma omp critical
    MPI_Get(local_trg, 1, dtype_trg, src_rank, disp, 1, dtype_src, win);

    // free the types
    MPI_Type_free(&dtype_trg);
    MPI_Type_free(&dtype_src);
    //-------------------------------------------------------------------------
}