#include "interpolator.hpp"

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
void Interpolator::Copy(const level_t dlvl, const lid_t shift[3], MemLayout *block_src, real_p data_src, MemLayout *block_trg, real_p data_trg) {
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
void Interpolator::Interpolate(const level_t dlvl, const lid_t shift[3], MemLayout *block_src, real_p data_src, MemLayout *block_trg, real_p data_trg) {
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
void Interpolator::Interpolate(const level_t dlvl, const lid_t shift[3], MemLayout *block_src, real_p data_src, MemLayout *block_trg, real_p data_trg, const real_t alpha, real_p data_cst) {
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
void Interpolator::DoMagic_(const level_t dlvl,const bool force_copy, const lid_t shift[3], MemLayout *block_src, real_p data_src, MemLayout *block_trg, real_p data_trg, const real_t alpha, real_p data_cst)
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
void Interpolator::Copy_(const level_t dlvl, const interp_ctx_t* ctx) {
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
void Interpolator::GetRma(const level_t dlvl, const lid_t shift[3], MemLayout *block_src, MPI_Aint disp_src, MemLayout *block_trg, real_t *data_trg, rank_t src_rank, MPI_Win win) {
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
    // get the corresponding MPI_Datatype for the source
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

void Interpolator::PutRma(const level_t dlvl, const lid_t shift[3], MemLayout *block_src, real_t *ptr_src, MemLayout *block_trg, MPI_Aint disp_trg, rank_t trg_rank, MPI_Win win) {
    m_assert(dlvl <= 1, "we cannot handle a difference in level > 1");
    m_assert(dlvl >= 0, "we cannot handle a level coarse ");
    m_assert(disp_trg >= 0, "the displacement is not positive: %ld", disp_trg);
    //-------------------------------------------------------------------------
    m_log("----------------- PUT RMA");
    m_log("entering interpolator with shift = %d %d %d", shift[0], shift[1], shift[2]);
    m_log("entering interpolator with srcstart = %d %d %d", block_src->start(0), block_src->start(1), block_src->start(2));
    m_log("entering interpolator with srcend = %d %d %d", block_src->end(0), block_src->end(1), block_src->end(2));
    m_log("entering interpolator with trgstart = %d %d %d", block_trg->start(0), block_trg->start(1), block_trg->start(2));
    m_log("entering interpolator with trgend = %d %d %d", block_trg->end(0), block_trg->end(1), block_trg->end(2));
    //................................................
    // get the corresponding MPI_Datatype for the target
    const lid_t  trg_start[3] = {block_trg->start(0), block_trg->start(1), block_trg->start(2)};
    const lid_t  trg_end[3]   = {block_trg->end(0), block_trg->end(1), block_trg->end(2)};
    MPI_Datatype dtype_trg    = ToMPIDatatype(trg_start, trg_end, block_trg->gs(), block_trg->stride(), 1);
    m_log("the trg = %d %d %d to %d %d %d", trg_start[0], trg_start[1], trg_start[2], trg_end[0], trg_end[1], trg_end[2]);

    //................................................
    // get the corresponding MPI_Datatype for the source
    const lid_t  scale        = (lid_t)pow(2, dlvl);
    const lid_t  src_start[3] = {shift[0] + block_trg->start(0) * scale, shift[1] + block_trg->start(1) * scale, shift[2] + block_trg->start(2) * scale};
    const lid_t  src_end[3]   = {shift[0] + block_trg->end(0) * scale, shift[1] + block_trg->end(1) * scale, shift[2] + block_trg->end(2) * scale};
    MPI_Datatype dtype_src    = ToMPIDatatype(src_start, src_end, block_src->gs(), block_src->stride(), scale);
    m_log("the src = %d %d %d to %d %d %d", src_start[0], src_start[1], src_start[2], src_end[0], src_end[1], src_end[2]);

    m_log("put to mirror num %d in rank %d",disp_trg/m_blockmemsize(1),trg_rank);

    //................................................
#ifndef NDEBUG
    int size_trg, size_src;
    MPI_Type_size(dtype_src, &size_src);
    MPI_Type_size(dtype_trg, &size_trg);
    // m_verb("src size = %d and the trg size = %d", size_src, size_trg);
    m_assert(size_trg == size_src, "the two sizes must match: src = %d vs trg = %d", size_src, size_trg);
#endif
    real_t * local_src = ptr_src + m_zeroidx(0, block_src) + m_midx(src_start[0], src_start[1], src_start[2], 0, block_src);
    MPI_Aint disp      = disp_trg + m_zeroidx(0, block_trg) + m_midx(trg_start[0], trg_start[1], trg_start[2], 0, block_trg);
#pragma omp critical
    // MPI_Get(local_trg, 1, dtype_trg, src_rank, disp, 1, dtype_src, win);
    MPI_Put(local_src, 1, dtype_src, trg_rank, disp, 1, dtype_trg, win);

    int size;
    MPI_Type_size(dtype_src, &size);
    m_log("wanna take %d bytes", size);
    real_t *tmp = reinterpret_cast<real_t *>(m_calloc(size));
    MPI_Status status;
    MPI_Sendrecv(local_src, 1, dtype_src, 0, 0, tmp, size / sizeof(real_t), M_MPI_REAL, 0, 0, MPI_COMM_SELF, &status);
    for (int iu = 0; iu < size / sizeof(real_t); iu++) {
        printf("%f ", tmp[iu]);
    }
    printf("\n ");
    m_free(tmp);

    // free the types
    MPI_Type_free(&dtype_trg);
    MPI_Type_free(&dtype_src);
    //-------------------------------------------------------------------------
}