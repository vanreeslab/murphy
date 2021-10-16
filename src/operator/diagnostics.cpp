#include "operator/diagnostics.hpp"

DetailVsError::DetailVsError(const Wavelet* interp) noexcept : BlockOperator(nullptr) {
    ghost_len_need_[0] = interp->nghost_front();
    ghost_len_need_[1] = interp->nghost_front();
}

void DetailVsError::operator()(const iter_t id, const std::string folder, const std::string suffix, const Grid* grid, Field* field, const lambda_error_t* sol) {
    m_begin;
    //--------------------------------------------------------------------------
    // allow to get only the detail distribution
    do_error_ = (sol != nullptr);

    //allocate the array
    size_t  size_blocks = (n_cat_ + 1) * (P8EST_MAXLEVEL) * ((do_error_) ? (n_cat_ + 1) : 1);
    n_blocks_loc_       = reinterpret_cast<bidx_t*>(m_calloc(size_blocks * sizeof(bidx_t)));
    n_blocks_glob_      = reinterpret_cast<bidx_t*>(m_calloc(size_blocks * sizeof(bidx_t)));

    // get the local distribution
    DoOpMesh(this, &DetailVsError::DoMagic, grid, grid->interp(), field, sol);

    // sum everything on rank 0 to dump
    const rank_t root = 0;
    m_assert(sizeof(bidx_t) == sizeof(int), "the two sizes must match");
    MPI_Reduce(n_blocks_loc_, n_blocks_glob_, size_blocks, MPI_INT, MPI_SUM, root, MPI_COMM_WORLD);

    // dump
    rank_t rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

#ifndef NDEBUG
    long nblock = grid->global_num_quadrants();
#endif

    FILE* file_diag;
    if (rank == root) {
        file_diag = fopen(std::string(folder + "/deterr_heat" + suffix + "_" + std::to_string(id) + ".data").c_str(), "a+");
#ifndef NDEBUG
        long block_count = 0;
#endif
        
        for (iter_t ierr = 0; ierr <= (do_error_*n_cat_); ++ierr) {
            for (level_t il=0; il<=P8EST_QMAXLEVEL;++il){
            for (iter_t idet = 0; idet <= n_cat_; ++idet) {
                // get the center of the category
                const real_t h_cat      = (log10(max_cat_) - log10(min_cat_)) / (n_cat_);
                const real_t err_bound = pow(10.0, log10(min_cat_) + ierr * h_cat);
                const real_t det_bound = pow(10.0, log10(min_cat_) + idet * h_cat);

                // write it down
                const size_t id = idet + (n_cat_ + 1) * (il + P8EST_MAXLEVEL * ierr);
                fprintf(file_diag, "%e;%d;%e;%d\n", err_bound * do_error_, il, det_bound, n_blocks_glob_[id]);
#ifndef NDEBUG
                block_count += n_blocks_glob_[id];
#endif
            }
            }
        }
        fclose(file_diag);
        m_assert(block_count == nblock, "the two counts must match: %ld %ld", block_count, nblock);
    }

    m_free(n_blocks_loc_);
    m_free(n_blocks_glob_);

    //--------------------------------------------------------------------------
    m_end;
}

void DetailVsError::DoMagic(const qid_t* qid, GridBlock* block, const Wavelet* interp, const Field* field, const lambda_error_t* sol) {
    m_assert(field->ghost_status(ghost_len_need_), "the ghost values of the velocity must be known! expected: %d %d, known: %d %d", ghost_len_need_[0], ghost_len_need_[1], field->get_ghost_len(0), field->get_ghost_len(1));
    //--------------------------------------------------------------------------
    // get the details for the block
    real_t        block_maxmin[2] = {0.0, 0.0};
    const MemSpan span_src(-ghost_len_need_[0], M_N + ghost_len_need_[1]);
    const MemSpan span_det = block->BlockSpan();
    ConstMemData  data_src = block->ConstData(field, 0);
    MemData       data_det(nullptr);
    interp->Details(&span_src, &data_src, &span_det, &data_det, 0.0, block_maxmin);

    // get the error for the block
    real_t erri = 0.0;

    auto op = [=, &erri](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
        // we need to discard the physical BC for the edges
        const real_t sol_val  = (*sol)(i0, i1, i2, block);
        const real_t data_val = data_src(i0, i1, i2);
        const real_t error    = data_val - sol_val;
        // update the block errors
        erri = m_max(std::fabs(error), erri);
    };
    if (do_error_) {
        for_loop(&op, span_det);
    }

    // get the category split
    const real_t h_cat = (log10(max_cat_) - log10(min_cat_)) / (n_cat_);

    iter_t  id_cat_det = m_max(0, m_min(n_cat_, (log10(block_maxmin[0]) - log10(min_cat_)) / h_cat));
    iter_t  id_cat_err = m_max(0, m_min(n_cat_, (log10(erri) - log10(min_cat_)) / h_cat));
    level_t id_level   = block->level();

    // got
    m_assert(n_cat_ >= id_cat_err && id_cat_err >= 0, "the cat id = %d must be >=0 and < %d", id_cat_err, n_cat_);
    m_assert(n_cat_ >= id_cat_det && id_cat_det >= 0, "the cat id = %d must be >=0 and < %d", id_cat_det, n_cat_);

    const size_t id = id_cat_det + (n_cat_ + 1) * (id_level + P8EST_MAXLEVEL * id_cat_err * do_error_);
    n_blocks_loc_[id] += 1;

    //--------------------------------------------------------------------------
};