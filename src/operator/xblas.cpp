#include "xblas.hpp"

#include "core/doop.hpp"
#include "core/forloop.hpp"

//------------------------------------------------------------------------------
BMax::BMax() noexcept : BlockOperator(nullptr){};
BMax::BMax(const bidx_t* ghost_len) noexcept : BlockOperator(ghost_len){};

/**
 * @brief computes the max(fabs(value)) over the grid
 */
real_t BMax::operator()(const ForestGrid* grid, const Field* fid_x) const {
    m_begin;
    m_assert(fid_x->ghost_status(ghost_len_need_), "The field  <%s> must have enough valid GP for the refinement - required %d %d, known %d %d", fid_x->name().c_str(), ghost_len_need_[0], ghost_len_need_[1], fid_x->get_ghost_len(0), fid_x->get_ghost_len(1));
    //--------------------------------------------------------------------------
    real_t max_grid = 0.0;
    // go on the blocks
    DoOpMesh(this, &BMax::ComputeBMaxGridBlock, grid, fid_x, &max_grid);
    // update the ghost - not needed
    // allreduce sync:
    real_t max_global = 0.0;
    MPI_Allreduce(&max_grid, &max_global, 1, M_MPI_REAL, MPI_MAX, MPI_COMM_WORLD);
    //--------------------------------------------------------------------------
    m_end;
    return max_global;
}

void BMax::ComputeBMaxGridBlock(const qid_t* qid, const CartBlock* block, const Field* fid_x, real_t* max) const {
    //--------------------------------------------------------------------------
    // get the data for each direction
    real_t local_max = 0.0;
    for (lda_t ida = 0; ida < fid_x->lda(); ++ida) {
        // for the current ida
        const ConstMemData data = block->data(fid_x, ida);

        auto op = [=, &local_max](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
            local_max = m_max(fabs(data(i0, i1, i2)), local_max);
        };
        for_loop(&op, span_);
    }
    max[0] = m_max(local_max, max[0]);
    //--------------------------------------------------------------------------
}

//------------------------------------------------------------------------------
BMinMax::BMinMax() noexcept : BlockOperator(nullptr){};
BMinMax::BMinMax(const bidx_t* ghost_len) noexcept : BlockOperator(ghost_len){};

/**
 * @brief compute the min and the max of fid_x dimension per dimension and store them in min and max (must have the dimension size)
 */
void BMinMax::operator()(const ForestGrid* grid, const Field* fid_x, real_t* min, real_t* max) const {
    m_begin;
    m_assert(fid_x->ghost_status(ghost_len_need_), "the field <%s> must be up to date", fid_x->name().c_str());
    //--------------------------------------------------------------------------
    // go on the blocks, for each dim separately
    for (lda_t ida = 0; ida < fid_x->lda(); ida++) {
        real_t res[2] = {std::numeric_limits<real_t>::max(),
                         std::numeric_limits<real_t>::min()};

        // store the dimension and go!
        DoOpMesh(this, &BMinMax::ComputeBMinMaxGridBlock, grid, fid_x, ida, res);

        // update the ghost - not needed
        // allreduce sync:
        MPI_Allreduce(res + 0, min + ida, 1, M_MPI_REAL, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(res + 1, max + ida, 1, M_MPI_REAL, MPI_MAX, MPI_COMM_WORLD);
    }
    //--------------------------------------------------------------------------
    m_end;
}

void BMinMax::ComputeBMinMaxGridBlock(const qid_t* qid, const CartBlock* block, const Field* fid_x, const lda_t ida, real_t res[2]) const {
    //-------------------------------------------------------------------------
    const ConstMemData data = block->data(fid_x, ida);
    // m_assume_aligned(data);

    auto op = [=, &res](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
        const real_t val = data(i0, i1, i2);
        // reduce the min/max
        res[0] = m_min(val, res[0]);
        res[1] = m_max(val, res[1]);
    };

    for_loop(&op, span_);
    //-------------------------------------------------------------------------
}
//==============================================================================
// //------------------------------------------------------------------------------
// BMoment::BMoment() noexcept : BlockOperator(nullptr) {
//     //--------------------------------------------------------------------------
//     // we need one more point at the back of the block
//     ghost_len_need_[0] = ghost_len_res_[0] + 2;  // + interp->nghost_front_refine();
//     ghost_len_need_[1] = ghost_len_res_[1] + 3;  // + interp->nghost_back_refine();

//     //--------------------------------------------------------------------------
// };
// BMoment::BMoment(const bidx_t* ghost_len) noexcept : BlockOperator(ghost_len) {
//     //--------------------------------------------------------------------------
//     // get the ghost sizes right
//     ghost_len_need_[0] = ghost_len_res_[0] + 2;  // interp->nghost_front_refine();
//     ghost_len_need_[1] = ghost_len_res_[1] + 3;  // interp->nghost_back_refine();
//     //--------------------------------------------------------------------------
// };

/**
 * @brief compute the 0th and the first moment of the block
 * 
 * @param moment0 the value of the 0th moment = integral (length = #lda of the field)
 * @param moment1 the vlaue of the 1st moments: m_x, m_y and m_z (length = 3 * #lda of the field)
 */
void BMoment::operator()(const ForestGrid* grid, const Field* fid_x, real_t* moment0, real_t* moment1) const {
    m_begin;
    m_assert(IsGhostValid(fid_x), "the field <%s> must be up to date", fid_x->name().c_str());
    //--------------------------------------------------------------------------
    for (lda_t ida = 0; ida < fid_x->lda(); ida++) {
        // #ifdef M_MPI_AGGRESSIVE
        //         MPI_Request requests[4];
        // #endif
        // moment 0
        auto op_moment_0 = [=](const bidx_t i0, const bidx_t i1, const bidx_t i2, const CartBlock* block) -> real_t {
            ConstMemData data_src = block->ConstData(fid_x, ida);
            return data_src(i0, i1, i2);
        };

        // #ifdef M_MPI_AGGRESSIVE
        //         ComputeIntegral(grid, op_moment_0, moment0 + ida, requests + 0);
        // #else
        ComputeIntegral(grid, op_moment_0, moment0 + ida);
        // #endif

        for (lda_t dir = 0; dir < 3; ++dir) {
            // moment 0
            auto op_moment_1 = [=](const bidx_t i0, const bidx_t i1, const bidx_t i2, const CartBlock* block) -> real_t {
                real_t pos[3];
                block->pos(i0, i1, i2, pos);
                ConstMemData data_src = block->ConstData(fid_x, ida);
                return data_src(i0, i1, i2) * pos[dir];
            };
            // #ifdef M_MPI_AGGRESSIVE
            //             ComputeIntegral(grid, op_moment_1, moment1+ dir + 3*ida, requests + 1 + dir);
            // #else
            ComputeIntegral(grid, op_moment_1, moment1 + dir + 3 * ida);
            // #endif
        }
        // #ifdef M_MPI_AGGRESSIVE
        //         MPI_Waitall(4, requests, MPI_STATUSES_IGNORE);
        // #endif
        m_log("moments are %e %e %e %e", moment0[0], moment1[0], moment1[1], moment1[2]);
    }
    //--------------------------------------------------------------------------
    m_end;
}

// void BMoment::ComputeBMomentGridBlock(const qid_t* qid, const CartBlock* block, const Field* fid_x, const lda_t ida, real_t moments[4]) const {
//     //--------------------------------------------------------------------------
//     const real_t* h        = block->hgrid();
//     ConstMemData  data_src = block->data(fid_x, ida);

//     real_t lmoment0    = 0.0;
//     real_t lmoment1[3] = {0.0, 0.0, 0.0};

//     // let's go!
//     auto op = [=, &lmoment0, &lmoment1](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
//         // simpson weights
//         const bidx_t n_loop =4 ;
//         const bidx_t shift[3] = {1,1,1};
//         const real_t weights[4] = {-1.0 / 24.0, 13.0 / 24.0, 13.0 / 24.0, -1.0 / 24.0};
//         // const bidx_t n_loop     = 2;
//         // const bidx_t shift[3]   = {0, 0, 0};
//         // const real_t weights[2] = {0.5, 0.5};
//         // local data
//         LocalData ldata(&data_src, i0, i1, i2);

//         // get the position
//         real_t origin[3];
//         block->pos(i0, i1, i2, origin);

//         real_t local_m00 = 0.0;
//         real_t local_m10 = 0.0;
//         real_t local_m11 = 0.0;
//         real_t local_m12 = 0.0;

// #ifndef NDEBUG
//         real_t sum_fact = 0.0;
// #endif
//         for (bidx_t is2 = 0; is2 < n_loop; ++is2) {
//             for (bidx_t is1 = 0; is1 < n_loop; ++is1) {
//                 for (bidx_t is0 = 0; is0 < n_loop; ++is0) {
//                     const real_t fact  = weights[is0] * weights[is1] * weights[is2];
//                     const real_t value = ldata(is0 - shift[0], is1 - shift[1], is2 - shift[2]);
//                     local_m00 += fact * value;
//                     local_m10 += fact * value * (origin[0] + h[0] * (is0 - shift[0]));
//                     local_m11 += fact * value * (origin[1] + h[1] * (is1 - shift[1]));
//                     local_m12 += fact * value * (origin[2] + h[2] * (is2 - shift[2]));

// #ifndef NDEBUG
//                     sum_fact += fact;
// #endif
//                 }
//             }
//         }
//         m_assert(std::fabs(sum_fact - 1.0) < 10.0 * std::numeric_limits<real_t>::epsilon(), "the sum must be 1 -> error = %e", std::fabs(sum_fact - 1.0));
//         lmoment0 += local_m00;
//         lmoment1[0] += local_m10;
//         lmoment1[1] += local_m11;
//         lmoment1[2] += local_m12;
//     };

//     m_assert(((span_.end[0] - span_.start[0]) % 2) == 0, "the span must be a modulo of 2: from %d to %d", span_.start[0], span_.end[0]);
//     m_assert(((span_.end[1] - span_.start[1]) % 2) == 0, "the span must be a modulo of 2: from %d to %d", span_.start[1], span_.end[1]);
//     m_assert(((span_.end[2] - span_.start[2]) % 2) == 0, "the span must be a modulo of 2: from %d to %d", span_.start[2], span_.end[2]);

//     // create a new span that will take the last GP
//     // const bidx_t span_shift[2][3] = {{0, 0, 0}, {-1, -1, -1}};
//     // MemSpan      moment_span(&span_, span_shift);
//     // const bidx_t jumps[3] = {1, 1, 1};
//     for_loop(&op, span_);

//     const real_t vol = block->hgrid(0) * block->hgrid(1) * block->hgrid(2);
//     moments[0] += vol * lmoment0;
//     moments[1] += vol * lmoment1[0];
//     moments[2] += vol * lmoment1[1];
//     moments[3] += vol * lmoment1[2];
//     //--------------------------------------------------------------------------
// };

//------------------------------------------------------------------------------
BAvg::BAvg() noexcept : BlockOperator(nullptr){};
BAvg::BAvg(const bidx_t* ghost_len) noexcept : BlockOperator(ghost_len){};

/**
 * @brief compute the sum of the values on the grid
 */
void BAvg::operator()(const ForestGrid* grid, const Field* fid_x, real_t* sum_global) const {
    m_begin;
    //--------------------------------------------------------------------------
    // reset the sum
    real_t sum_local = 0.0;
    // go on the blocks, for each dim separately
    for (lda_t ida = 0; ida < fid_x->lda(); ++ida) {
        DoOpMesh(this, &BAvg::ComputeBAvgGridBlock, grid, fid_x, ida, &sum_local);
    }
    // allreduce sync:
    MPI_Allreduce(&sum_local, sum_global, 1, M_MPI_REAL, MPI_SUM, MPI_COMM_WORLD);
    //--------------------------------------------------------------------------
    m_end;
}

void BAvg::ComputeBAvgGridBlock(const qid_t* qid, const CartBlock* block, const Field* fid_x, const lda_t ida, real_t* sum) const {
    //--------------------------------------------------------------------------
    // get the starting pointer:
    const real_t*      h    = block->hgrid();
    const ConstMemData data = block->data(fid_x, ida);

    real_t sum_local = 0.0;
    // let's go!
    auto op = [=, &sum_local](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
        sum_local += data(i0, i1, i2);
    };
    for_loop(&op, span_);
    sum[0] += sum_local * h[0] * h[1] * h[2];
    //--------------------------------------------------------------------------
}

//------------------------------------------------------------------------------
BDensity::BDensity() noexcept : BlockOperator(nullptr){};
BDensity::BDensity(const bidx_t* ghost_len) noexcept : BlockOperator(ghost_len){};

/**
 * @brief compute the sum of the values on the grid
 */
void BDensity::operator()(const ForestGrid* grid, real_t* density_global) const {
    m_begin;
    //--------------------------------------------------------------------------
    // reset the sum
    real_t density_local = 0.0;
    // go on the blocks, for each dim separately
    DoOpMesh(this, &BDensity::ComputeBDensityGridBlock, grid, &density_local);
    // allreduce sync:
    MPI_Allreduce(&density_local, density_global, 1, M_MPI_REAL, MPI_SUM, MPI_COMM_WORLD);
    density_global[0] = sqrt(density_global[0]);
    //--------------------------------------------------------------------------
    m_end;
}

void BDensity::ComputeBDensityGridBlock(const qid_t* qid, const CartBlock* block, real_t* density) const {
    //--------------------------------------------------------------------------
    // get the starting pointer:
    const real_t level  = block->level();
    const real_t length = block->length();
    density[0] += (level * level) * (length * length * length);
    //--------------------------------------------------------------------------
}

//------------------------------------------------------------------------------
// BDiscreteMoment::BDiscreteMoment() noexcept : BlockOperator(nullptr){};
// BDiscreteMoment::BDiscreteMoment(const Wavelet*  interp) noexcept : BlockOperator(interp){};

// /**
//  * @brief compute the 0th and the first moment of the block
//  *
//  * @param moment0 the value of the 0th moment = integral (length = #lda of the field)
//  * @param moment1 the vlaue of the 1st moments: m_x, m_y and m_z (length = 3 * #lda of the field)
//  */
// void BDiscreteMoment::operator()(const ForestGrid*  grid, const Field*  fid_x, real_t* moment0, real_t* moment1) {
//     m_begin;
//     //--------------------------------------------------------------------------
//     // go on the blocks, for each dim separately
//     for (lda_t ida = 0; ida < fid_x->lda(); ida++) {
//         // reset the moments
//         moment0_    = 0.0;  // moment0 = integral
//         moment1_[0] = 0.0;  // moment_x
//         moment1_[1] = 0.0;  // moment_y
//         moment1_[2] = 0.0;  // moment_z

//         // store the dimension and go!
//         ida_ = ida;
//         DoOpMesh(this, &BDiscreteMoment::ComputeBDiscreteMomentGridBlock, grid, fid_x);

//         // update the ghost - not needed
//         // allreduce sync:
//         MPI_Allreduce(&moment0_, moment0 + ida, 1, M_MPI_REAL, MPI_SUM, MPI_COMM_WORLD);
//         MPI_Allreduce(&moment1_, moment1 + 3 * ida, 3, M_MPI_REAL, MPI_SUM, MPI_COMM_WORLD);
//     }
//     //--------------------------------------------------------------------------
//     m_end;
// }

// /**
//  * @brief Integrate the discrete moments on the block
//  *
//  * @param qid
//  * @param block
//  * @param fid_x
//  */
// void BDiscreteMoment::ComputeBDiscreteMomentGridBlock(const qid_t*  qid, GridBlock*  block, const Field*  fid_x) {
//     m_assert((end_ - start_) % 4 == 0, "the span done = %d to %d must be a modulo of 4", start_, end_);
//     //--------------------------------------------------------------------------
//     // get the starting pointer:
//     const real_t* data = block->data(fid_x, ida_).Read();

//     real_t lmoment0    = 0.0;
//     real_t lmoment1[3] = {0.0, 0.0, 0.0};

//     // let's go!
//     for (bidx_t i2 = start_; i2 < end_; ++i2) {
//         for (bidx_t i1 = start_; i1 < end_; ++i1) {
//             for (bidx_t i0 = start_; i0 < end_; ++i0) {
//                 const real_t f = data[m_idx(i0, i1, i2)];

//                 // get the position
//                 real_t pos[3];
//                 m_pos(pos, i0, i1, i2, block->hgrid(), block->xyz());

//                 // get the coefficient

//                 lmoment0 += f;
//                 lmoment1[0] += f * pos[0];
//                 lmoment1[1] += f * pos[1];
//                 lmoment1[2] += f * pos[2];
//             }
//         }
//     }

//     const real_t vol = block->hgrid(0) * block->hgrid(1) * block->hgrid(2);
//     moment0_ += vol * lmoment0;
//     moment1_[0] += vol * lmoment1[0];
//     moment1_[1] += vol * lmoment1[1];
//     moment1_[2] += vol * lmoment1[2];
//     //--------------------------------------------------------------------------
// }
