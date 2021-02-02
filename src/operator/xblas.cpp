#include "xblas.hpp"
#include "core/doop.hpp"

//-----------------------------------------------------------------------------
BMax::BMax()noexcept : BlockOperator(nullptr){};
BMax::BMax(m_ptr<const Wavelet> interp)noexcept : BlockOperator(interp){};

real_t BMax::operator()(m_ptr<const ForestGrid> grid, m_ptr<const Field> fid_x) {
    m_begin;
    //-------------------------------------------------------------------------
    max_ = std::numeric_limits<real_t>::min();
    // go on the blocks
    DoOpMesh(this, &BMax::ComputeBMaxGridBlock, grid, fid_x);
    // update the ghost - not needed
    // allreduce sync:
    real_t max_global = 0.0;
    MPI_Allreduce(&max_, &max_global, 1, M_MPI_REAL, MPI_MAX, MPI_COMM_WORLD);
    //-------------------------------------------------------------------------
    m_end;
    return max_global;
}

void BMax::ComputeBMaxGridBlock(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, m_ptr<const Field> fid_x) {
    //-------------------------------------------------------------------------
    const sid_t lda = fid_x->lda();
    for (sid_t ida = 0; ida < lda; ida++) {
        // get the data pointers
        const real_t* data_x = block->data(fid_x, ida).Read();
        m_assume_aligned(data_x);
        // get the correct place given the current thread and the dimension
        for (bidx_t i2 = start_; i2 < end_; i2++) {
            for (bidx_t i1 = start_; i1 < end_; i1++) {
                for (bidx_t i0 = start_; i0 < end_; i0++) {
                    const size_t idx = m_idx(i0, i1, i2);
                    max_             = m_max(fabs(data_x[idx]), max_);
                }
            }
        }
    }
    //-------------------------------------------------------------------------
}


//-----------------------------------------------------------------------------
BMinMax::BMinMax()noexcept : BlockOperator(nullptr){};
BMinMax::BMinMax(m_ptr<const Wavelet> interp)noexcept : BlockOperator(interp){};

/**
 * @brief compute the min and the max of fid_x dimension per dimension and store them in min and max (must have the dimension size)
 */
void BMinMax::operator()(m_ptr<const ForestGrid> grid, m_ptr<const Field> fid_x, real_t* min, real_t* max) {
    m_begin;
    //-------------------------------------------------------------------------
    // go on the blocks, for each dim separately
    for (lda_t ida = 0; ida < fid_x->lda(); ida++) {
        max_ = std::numeric_limits<real_t>::min();
        min_ = std::numeric_limits<real_t>::max();

        // store the dimension and go!
        ida_ = ida;
        DoOpMesh(this, &BMinMax::ComputeBMinMaxGridBlock, grid, fid_x);

        // update the ghost - not needed
        // allreduce sync:
        MPI_Allreduce(&min_, min + ida, 1, M_MPI_REAL, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&max_, max + ida, 1, M_MPI_REAL, MPI_MAX, MPI_COMM_WORLD);
    }
    //-------------------------------------------------------------------------
    m_end;
}

void BMinMax::ComputeBMinMaxGridBlock(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, m_ptr<const Field> fid_x) {
    //-------------------------------------------------------------------------
    const real_t* data_x = block->data(fid_x, ida_).Read();
    m_assume_aligned(data_x);

    for (bidx_t i2 = start_; i2 < end_; i2++) {
        for (bidx_t i1 = start_; i1 < end_; i1++) {
            for (bidx_t i0 = start_; i0 < end_; i0++) {
                const size_t idx = m_idx(i0, i1, i2);
                min_             = m_min(data_x[idx], min_);
                max_             = m_max(data_x[idx], max_);
            }
        }
    }
    //-------------------------------------------------------------------------
}

//-----------------------------------------------------------------------------
BMoment::BMoment() noexcept: BlockOperator(nullptr){};
BMoment::BMoment(m_ptr<const Wavelet> interp) noexcept: BlockOperator(interp){};

/**
 * @brief compute the 0th and the first moment of the block
 * 
 * @param moment0 the value of the 0th moment = integral (length = #lda of the field)
 * @param moment1 the vlaue of the 1st moments: m_x, m_y and m_z (length = 3 * #lda of the field)
 */
void BMoment::operator()(m_ptr<const ForestGrid> grid, m_ptr<const Field> fid_x, real_t* moment0, real_t* moment1) {
    m_begin;
    m_assert(fid_x->ghost_status(), "the field <%s> must be up to date", fid_x->name().c_str());
    //-------------------------------------------------------------------------
    // go on the blocks, for each dim separately
    for (lda_t ida = 0; ida < fid_x->lda(); ida++) {
        // reset the moments
        moment0_    = 0.0;  // moment0 = integral
        moment1_[0] = 0.0;  // moment_x
        moment1_[1] = 0.0;  // moment_y
        moment1_[2] = 0.0;  // moment_z

        // store the dimension and go!
        ida_ = ida;
        DoOpMesh(this, &BMoment::ComputeBMomentGridBlock, grid, fid_x);

        // update the ghost - not needed
        // allreduce sync:
        MPI_Allreduce(&moment0_, moment0 + ida, 1, M_MPI_REAL, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&moment1_, moment1 + 3 * ida, 3, M_MPI_REAL, MPI_SUM, MPI_COMM_WORLD);
    }
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief Integrate the different moments on the block, using Simpson 3/8 in 3D
 * 
 * cfr. Abramowitz p886, formula 25.4.14:
 *      int_x0^x4 f(x) dx = 2/45 * h * (7 f0 + 32 f1 + 12 f2 + 32 f3 + 7 f4)
 * 
 * @param qid 
 * @param block 
 * @param fid_x 
 */
void BMoment::ComputeBMomentGridBlock(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, m_ptr<const Field> fid_x) {
    m_assert((end_ - start_) % 4 == 0, "the span done = %d to %d must be a modulo of 4", start_, end_);
    //-------------------------------------------------------------------------
    // get the starting pointer:
    const real_t* data = block->data(fid_x, ida_).Read();

    real_t lmoment0    = 0.0;
    real_t lmoment1[3] = {0.0, 0.0, 0.0};

    // let's go!
    for (bidx_t i2 = start_; i2 < end_; i2 += 4) {
        for (bidx_t i1 = start_; i1 < end_; i1 += 4) {
            for (bidx_t i0 = start_; i0 < end_; i0 += 4) {
                const real_t* f = data + m_idx(i0, i1, i2);

                // define the simpson coefficients - cfr Abramowitz p887, formula 25.4.13
                constexpr real_t c[5] = {7.0, 32.0, 12.0, 32.0, 7.0};

                // integrate simpson in 3D
                for (bidx_t s2 = 0; s2 <= 4; ++s2) {
                    for (bidx_t s1 = 0; s1 <= 4; ++s1) {
                        for (bidx_t s0 = 0; s0 <= 4; ++s0) {
                            // get the local value
                            const real_t value = f[m_idx(s0, s1, s2)];

                            // get the position
                            real_t pos[3];
                            m_pos(pos, i0 + s0, i1 + s1, i2 + s2, block->hgrid(), block->xyz());

                            // get the coefficient
                            const real_t coef = (c[s0] * c[s1] * c[s2]);

                            lmoment0 += coef * value;
                            lmoment1[0] += coef * value * pos[0];
                            lmoment1[1] += coef * value * pos[1];
                            lmoment1[2] += coef * value * pos[2];
                        }
                    }
                }
            }
        }
    }
    const real_t vol = block->hgrid(0) * block->hgrid(2) * block->hgrid(2);
    moment0_ += pow(2.0 / 45.0, 3) * vol * lmoment0;
    moment1_[0] += pow(2.0 / 45.0, 3) * vol * lmoment1[0];
    moment1_[1] += pow(2.0 / 45.0, 3) * vol * lmoment1[1];
    moment1_[2] += pow(2.0 / 45.0, 3) * vol * lmoment1[2];
    //-------------------------------------------------------------------------
}

//-----------------------------------------------------------------------------
BDiscreteMoment::BDiscreteMoment() noexcept : BlockOperator(nullptr){};
BDiscreteMoment::BDiscreteMoment(m_ptr<const Wavelet> interp) noexcept : BlockOperator(interp){};

/**
 * @brief compute the 0th and the first moment of the block
 * 
 * @param moment0 the value of the 0th moment = integral (length = #lda of the field)
 * @param moment1 the vlaue of the 1st moments: m_x, m_y and m_z (length = 3 * #lda of the field)
 */
void BDiscreteMoment::operator()(m_ptr<const ForestGrid> grid, m_ptr<const Field> fid_x, real_t* moment0, real_t* moment1) {
    m_begin;
    //-------------------------------------------------------------------------
    // go on the blocks, for each dim separately
    for (lda_t ida = 0; ida < fid_x->lda(); ida++) {
        // reset the moments
        moment0_    = 0.0;  // moment0 = integral
        moment1_[0] = 0.0;  // moment_x
        moment1_[1] = 0.0;  // moment_y
        moment1_[2] = 0.0;  // moment_z

        // store the dimension and go!
        ida_ = ida;
        DoOpMesh(this, &BDiscreteMoment::ComputeBDiscreteMomentGridBlock, grid, fid_x);

        // update the ghost - not needed
        // allreduce sync:
        MPI_Allreduce(&moment0_, moment0 + ida, 1, M_MPI_REAL, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&moment1_, moment1 + 3 * ida, 3, M_MPI_REAL, MPI_SUM, MPI_COMM_WORLD);
    }
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief Integrate the different moments on the block, using Simpson 3/8 in 3D
 * 
 * cfr. Abramowitz p886, formula 25.4.14:
 *      int_x0^x4 f(x) dx = 2/45 * h * (7 f0 + 32 f1 + 12 f2 + 32 f3 + 7 f4)
 * 
 * @param qid 
 * @param block 
 * @param fid_x 
 */
void BDiscreteMoment::ComputeBDiscreteMomentGridBlock(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, m_ptr<const Field> fid_x) {
    m_assert((end_ - start_) % 4 == 0, "the span done = %d to %d must be a modulo of 4", start_, end_);
    //-------------------------------------------------------------------------
    // get the starting pointer:
    const real_t* data = block->data(fid_x, ida_).Read();

    real_t lmoment0    = 0.0;
    real_t lmoment1[3] = {0.0, 0.0, 0.0};

    // let's go!
    for (bidx_t i2 = start_; i2 < end_; ++i2) {
        for (bidx_t i1 = start_; i1 < end_; ++i1) {
            for (bidx_t i0 = start_; i0 < end_; ++i0) {
                const real_t f = data[m_idx(i0, i1, i2)];

                // get the position
                real_t pos[3];
                m_pos(pos, i0, i1, i2, block->hgrid(), block->xyz());

                // get the coefficient

                lmoment0 += f;
                lmoment1[0] += f * pos[0];
                lmoment1[1] += f * pos[1];
                lmoment1[2] += f * pos[2];
            }
        }
    }

    const real_t vol = block->hgrid(0) * block->hgrid(2) * block->hgrid(2);
    moment0_ += vol * lmoment0;
    moment1_[0] += vol * lmoment1[0];
    moment1_[1] += vol * lmoment1[1];
    moment1_[2] += vol * lmoment1[2];
    //-------------------------------------------------------------------------
}
