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
static inline MPI_Datatype ToMPIDatatype_ptr(const lid_t start[3], const lid_t end[3], const lid_t gs, const lid_t stride, const lid_t scale) {
        m_begin;
        //-------------------------------------------------------------------------
        // convert the
        lid_t start_g[3] = {start[0] + gs, start[1] + gs, start[2] + gs};
        lid_t end_g[3]   = {end[0] + gs, end[1] + gs, end[2] + gs};

        // MPI_Aint     lb, extent;
        MPI_Datatype x_type, xy_type, xyz_type;
        // MPI_Datatype x_tmp, xy_tmp, xyz_tmp;

        //................................................
        // do x type
        int      count_x  = end_g[0] - start_g[0];
        MPI_Aint stride_x = sizeof(real_t) * scale;
        m_assert(count_x > 0, "we at least need to take 1 element");
        m_assert(count_x <= stride, "we cannot take more element than the stride");
        // MPI_Type_create_hvector(count_x, 1, stride_x, M_MPI_REAL, &x_tmp);
        MPI_Type_create_hvector(count_x, 1, stride_x, M_MPI_REAL, &x_type);
        // MPI_Type_commit(&x_type);

        // shift to reach the correct spot
        // MPI_Type_get_extent(x_tmp, &lb, &extent);
        // lb = MPI_Aint_add(lb, stride_x * start_g[0]);
        // MPI_Type_create_resized(x_tmp, lb, extent, &x_type);
        // MPI_Type_commit(&x_type);
        // m_log("x type: bound = %ld, extend = %ld, count = %d and stride = %d",lb,extent,count_x,stride_x);

        //................................................
        // do y type
        int      count_y  = end_g[1] - start_g[1];
        MPI_Aint stride_y = sizeof(real_t) * stride * scale;
        m_assert(count_y > 0, "we at least need to take 1 element");
        m_assert(count_y <= stride, "we cannot take more element than the stride");
        // MPI_Type_create_hvector(count_y, 1, stride_y, x_type, &xy_tmp);
        // MPI_Type_commit(&xy_tmp);
        MPI_Type_create_hvector(count_y, 1, stride_y, x_type, &xy_type);
        // MPI_Type_commit(&xy_type);

        // shift to reach the correct spot
        // MPI_Type_get_extent(xy_tmp, &lb, &extent);
        // lb = MPI_Aint_add(lb, stride_y * start_g[1]);
        // MPI_Type_create_resized(xy_tmp, lb, extent, &xy_type);
        // MPI_Type_commit(&xy_type);

        // m_log("xy type: bound = %ld, extend = %ld, count = %d and stride = %d",lb,extent,count_y,stride_y);

        //................................................
        // do z type
        int      count_z  = end_g[2] - start_g[2];
        MPI_Aint stride_z = sizeof(real_t) * stride * stride * scale;
        m_assert(count_z > 0, "we at least need to take 1 element");
        m_assert(count_z <= stride, "we cannot take more element than the stride");
        // MPI_Type_create_hvector(count_z, 1, stride_z, xy_type, &xyz_tmp);
        // MPI_Type_commit(&xyz_tmp);
        MPI_Type_create_hvector(count_z, 1, stride_z, xy_type, &xyz_type);
        MPI_Type_commit(&xyz_type);

        // shift to reach the correct spot
        // MPI_Type_get_extent(xyz_tmp, &lb, &extent);
        // lb = MPI_Aint_add(lb, stride_z * start_g[2]);
        // MPI_Type_create_resized(xyz_tmp, lb, extent, &xyz_type);
        // MPI_Type_commit(&xyz_type);

        // m_log("xyy type: bound = %ld, extend = %ld, count = %d and stride = %d",lb,extent,count_z,stride_z);

        //................................................
        // free the useless types and return
        MPI_Type_free(&x_type);
        MPI_Type_free(&xy_type);
        // MPI_Type_free(&x_tmp);
        // MPI_Type_free(&xy_tmp);
        // MPI_Type_free(&xyz_tmp);

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
void Interpolator::GetRma(const sid_t dlvl, const lid_t shift[3], MemLayout *block_src, MPI_Aint disp_src, MemLayout *block_trg, real_t* data_trg, int src_rank, MPI_Win win) {
    m_assert(dlvl <= 1, "we cannot handle a difference in level > 1");
    m_assert(dlvl >= 0, "we cannot handle a level coarse ");
    //-------------------------------------------------------------------------
    m_log("-----------------");
    m_log("entering interpolator with shift = %d %d %d", shift[0], shift[1], shift[2]);
    m_log("entering interpolator with srcstart = %d %d %d", block_src->start(0), block_src->start(1), block_src->start(2));
    m_log("entering interpolator with srcend = %d %d %d", block_src->end(0), block_src->end(1), block_src->end(2));
    m_log("entering interpolator with trgstart = %d %d %d", block_trg->start(0), block_trg->start(1), block_trg->start(2));
    m_log("entering interpolator with trgend = %d %d %d", block_trg->end(0), block_trg->end(1), block_trg->end(2));

    //................................................
    // get the corresponding MPI_Datatype for the target
    const lid_t  trg_start[3] = {block_trg->start(0), block_trg->start(1), block_trg->start(2)};
    const lid_t  trg_end[3]   = {block_trg->end(0), block_trg->end(1), block_trg->end(2)};
    MPI_Datatype dtype_trg    = ToMPIDatatype_ptr(trg_start, trg_end, block_trg->gs(), block_trg->stride(), 1);

    m_log("the trg = %d %d %d to %d %d %d",trg_start[0],trg_start[1],trg_start[2],trg_end[0],trg_end[1],trg_end[2]);


    //................................................
    // get the corresponding MPI_Datatype for the target
    const lid_t  scale        = (lid_t)pow(2, dlvl);
    const lid_t  src_start[3] = {shift[0] + block_trg->start(0) / scale, shift[1] + block_trg->start(1) / scale, shift[2] + block_trg->start(2) / scale};
    const lid_t  src_end[3]   = {shift[0] + block_trg->end(0) / scale, shift[1] + block_trg->end(1) / scale, shift[2] + block_trg->end(2) / scale};
    MPI_Datatype dtype_src    = ToMPIDatatype_ptr(src_start, src_end, block_src->gs(), block_src->stride(), scale);

    m_log("the src = %d %d %d to %d %d %d", src_start[0], src_start[1], src_start[2], src_end[0], src_end[1], src_end[2]);

    // int size;
    // MPI_Type_size(dtype_trg, &size);

    // // m_log("trg_size = %d", size);
    // // real_t *temp = reinterpret_cast<real_t *>(m_calloc(size));
    // // size         = size / sizeof(real_t);
    // // MPI_Status status;
    // // MPI_Sendrecv(ptr_trg, 1, dtype_trg, 0, 0, temp, size, M_MPI_REAL, 0, 0, MPI_COMM_SELF, &status);
    // for (int ib = 0; ib < size; ib++) {
    //     temp[ib] = 1.0;
    //     printf("[%d]: %f ", ib, temp[ib]);
    // }
    // printf("\n");

    // int rank;
    // MPI_Comm_rank(MPI_COMM_SELF, &rank);
    // real_t * recv_ptr = ptr_trg + m_zeroidx(0, block_trg) + m_midx(trg_start[0], trg_start[1], trg_start[2], 0, block_trg);
    // MPI_Sendrecv(temp, size, M_MPI_REAL, rank, 0, recv_ptr, 1, dtype_trg, rank, 0, MPI_COMM_SELF, &status);
    // m_free(temp);

    // m_log("check");
    // for (int i2 = trg_start[2]; i2 < trg_end[2]; i2++) {
    //     for (int i1 = trg_start[1]; i1 < trg_end[1]; i1++) {
    //         for (int i0 = trg_start[0]; i0 < trg_end[0]; i0++) {
    //             printf("%f ", *(ptr_trg + m_zeroidx(0, block_trg) + m_midx(i0, i1, i2, 0, block_trg)));
    //         }
    //         printf("\n");
    //     }
    // }
    // MPI_Aint lb, extent;
    // MPI_Type_get_extent(dtype_trg,&lb,&extent);
    // real_t* point = ptr_trg+lb/sizeof(real_t);
    // printf("the target is")
    // for (lid_t ib=0; ib< extent/sizeof(real_t); ib++){
    //     printf
    // }

    // MPI_Aint lb, extent;
    // MPI_Type_get_true_extent(dtype_trg, &lb, &extent);
    // m_log("the true lb and extend of data_trg type is given by %ld %ld", lb, extent);
    // // MPI_Aint lb, extent;
    // MPI_Type_get_extent(dtype_trg, &lb, &extent);
    // m_log("vs the normal one %ld %ld", lb, extent);

    // m_log("accessing the mirrors with disp = %ld = %ld mirror", disp_src, disp_src / m_blockmemsize(1));
    m_assert(disp_src >= 0, "the displacement is not positive: %ld", disp_src);
    real_t *local_trg = data_trg + m_midx(trg_start[0], trg_start[1], trg_start[2], 0, block_trg);
    // move the displacement to the correct position
    MPI_Aint disp = disp_src + m_zeroidx(0, block_src) + m_midx(src_start[0], src_start[1], src_start[2], 0, block_src);
#pragma omp critical
    MPI_Get(local_trg, 1, dtype_trg, src_rank, disp, 1, dtype_src, win);

    // free the types
    MPI_Type_free(&dtype_trg);
    MPI_Type_free(&dtype_src);

    // int num_int;
    // int num_address;
    // int num_datatype;
    // int combiner;
    // MPI_Type_get_envelope(dtype_trg,&num_int,&num_address,&num_datatype,&combiner);
    //     // m_log("the trg datatype: num_in = %d, num_address = %d, num_datatype = %d, combiner = %d",num_int,num_address,num_datatype,combiner);

    //     lid_t start_g[3] = {block_trg->start(0) + block_src->gs(), block_trg->start(1) + block_src->gs(), block_trg->start(2) + block_src->gs()};
    //     lid_t end_g[3]   = {block_trg->end(0) + block_src->gs(), block_trg->end(1) + block_src->gs(), block_trg->end(2) + block_src->gs()};

    //     //................................................
    //     // the source datatype is driven by the number of dlvl
    //     lid_t scale = (lid_t)pow(2, dlvl);

    //     // create the datatype of one line in x -> a "block"
    //     MPI_Datatype xtype;
    //     int          count = (end_g[0] - start_g[0]);
    //     m_log("we take %d points (scale = %d) ", count, scale);
    //     MPI_Type_create_hvector(count, 1, scale * sizeof(real_t), M_MPI_REAL, &xtype);
    //     MPI_Type_commit(&xtype);

    //     // // do some stupid checks
    //     // MPI_Aint x_lb, x_extent;
    //     // MPI_Type_get_extent(xtype, &x_lb, &x_extent);
    //     // int size1 = 0;
    //     // MPI_Type_size(xtype, &size1);
    //     // m_log("xtype lower bound = %ld and extend = %ld -> size = %d", x_lb, x_extent, size1);

    //     // put the blocks together into a single datatype
    //     MPI_Datatype dtype_tmp;
    //     count                = (end_g[1] - start_g[1]) * (end_g[2] - start_g[2]);
    //     MPI_Aint stride_byte = block_src->stride() * scale * sizeof(real_t);
    //     MPI_Type_create_hvector(count, 1, stride_byte, xtype, &dtype_tmp);
    //     m_log("we repeat the blocks %d times with a stride of %d doubles = %d bytes", count, block_src->stride() * scale, stride_byte);
    //     MPI_Type_commit(&dtype_tmp);
    //     MPI_Type_free(&xtype);

    //     //................................................
    //     // finally add the shift from the (0,0,0) location
    //     MPI_Datatype dtype_src;
    //     MPI_Aint     lb, extent;
    //     MPI_Type_get_extent(dtype_tmp, &lb, &extent);
    //     m_log("mpi type_tmp lower bound = %ld and extend = %ld", lb, extent);
    //     // update the lower bound
    //     // lb = MPI_Aint_add(lb, (block_src->start(0) / scale) * (block_src->start(1) / scale) * (block_src->start(2) / scale) * sizeof(real_t));
    //     MPI_Aint offset_byte = sizeof(real_t) * m_sidx(shift[0] + start_g[0] / scale,
    //                                                    shift[1] + start_g[1] / scale,
    //                                                    shift[2] + start_g[2] / scale, 0, block_src->stride());
    //     lb                   = MPI_Aint_add(lb, offset_byte);
    //     m_log("shifting to byte number %ld",offset_byte);
    //     // m_log("shifting the type by %ld bytes = %d %d %d -> now lb = %d", m_sidx(start[0], start[1], start[2], 0, block_src->stride()), start[0], start[1], start[2], lb);
    //     // change the type and create a new one
    //     MPI_Type_create_resized(dtype_tmp, lb, extent, &dtype_src);
    //     MPI_Type_commit(&dtype_src);
    //     MPI_Type_free(&dtype_tmp);

    //     m_log("starting the get from proc %d to me: mirror number %ld", src_rank, disp_src / m_blockmemsize(1));
    //     // m_log("disp = %ld", disp_src);

    // // start the get in RMA language: origin = trg for us and target = src for us
    // #ifndef NDEBUG
    //     int src_size, trg_size;
    //     MPI_Type_size(dtype_src, &src_size);
    //     MPI_Type_size(dtype_trg, &trg_size);
    //     m_assert(src_size == trg_size, "the number of bytes should match in both the src and the trg");
    //         m_log("the src datatype = %d bytes and the trg = %d bytes", src_size, trg_size);
    // // #endif
    // // #pragma omp critical
    //     MPI_Get(ptr_trg, 1, dtype_trg, src_rank, disp_src, 1, dtype_src, win);

    //     // free the types
    //     MPI_Type_free(&dtype_trg);
    //     MPI_Type_free(&dtype_src);
    //-------------------------------------------------------------------------
}