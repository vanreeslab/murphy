#ifndef SRC_LAPLACIAN_IPP_
#define SRC_LAPLACIAN_IPP_

#include "laplacian.hpp"

/**
 * @brief computes the inner laplacian contribution of the dimension Stencil::ida_ of the source field on the traget field
 * 
 * @tparam length
 * @param qid
 * @param block 
 * @param fid_src the source field
 * @param fid_trg the target field
 */
template <sid_t length>
void LaplacianCross<length>::ApplyOpDerivInner(const qid_t* qid, GridBlock* block, const Field* fid_src, Field* fid_trg) {
    //-------------------------------------------------------------------------
    real_t* coef     = coef_ + (length / 2);
    real_t  scale[3] = {1.0/pow(block->hgrid(0), 2), 1.0/pow(block->hgrid(1), 2),1.0/ pow(block->hgrid(2), 2)};

    real_p data_src = block->data(fid_src, ida_);
    real_p data_trg = block->data(fid_trg, ida_);
    m_assume_aligned(data_trg);
    m_assume_aligned(data_src);

    // m_log("the coefs are %f %f %f %f %f",coef[-2],coef[-1],coef[0],coef[1],coef[2]);

    for (lid_t i2 = M_GS; i2 < M_N - M_GS; i2++) {
        for (lid_t i1 = M_GS; i1 < M_N - M_GS; i1++) {
            for (lid_t i0 = M_GS; i0 < M_N - M_GS; i0++) {
                // get the data pointer in front of the row for every cache line
                real_p       trg = data_trg + m_idx(i0, i1, i2);  // cache line for writting
                const real_p src = data_src + m_idx(i0, i1, i2);  // cache line for reading
                // loop on the stencil block
                trg[0] = 0.0;
                for (int is = -(length / 2); is <= (length / 2); is++) {
                    trg[0] += coef[is] * scale[0] * src[m_idx(is, 0, 0)];
                    trg[0] += coef[is] * scale[1] * src[m_idx(0, is, 0)];
                    trg[0] += coef[is] * scale[2] * src[m_idx(0, 0, is)];
                }
            }
        }
    }
    //-------------------------------------------------------------------------
}

/**
 * @brief computes the 6 different outer laplacian contribution of the dimension Stencil::ida_ of the source field's ghost values on the traget field
 * 
 * @tparam length 
 * @param qid 
 * @param block 
 * @param fid_src 
 * @param fid_trg 
 */
template <sid_t length>
void LaplacianCross<length>::ApplyOpDerivOuter(const qid_t* qid, GridBlock* block, const Field* fid_src, Field* fid_trg) {
    //-------------------------------------------------------------------------
    real_t* coef     = coef_ + (length / 2);
    real_t  scale[3] = {1.0/pow(block->hgrid(0), 2), 1.0/pow(block->hgrid(1), 2),1.0/ pow(block->hgrid(2), 2)};

    real_p data_src = block->data(fid_src, ida_);
    real_p data_trg = block->data(fid_trg, ida_);
    m_assume_aligned(data_trg);
    m_assume_aligned(data_src);

    //------ X -
    for (lid_t i2 = 0; i2 < M_N ; i2++) {
        for (lid_t i1 = 0; i1 < M_N; i1++) {
            for (lid_t i0 = 0; i0 < M_GS; i0++) {
                // get the data pointer in front of the row for every cache line
                real_p       trg = data_trg + m_idx(i0, i1, i2);  // cache line for writting
                const real_p src = data_src + m_idx(i0, i1, i2);  // cache line for reading
                // loop on the stencil block
                trg[0] = 0.0;
                for (int is = -(length / 2); is <= (length / 2); is++) {
                    trg[0] += coef[is] * scale[0] * src[m_idx(is, 0, 0)];
                    trg[0] += coef[is] * scale[1] * src[m_idx(0, is, 0)];
                    trg[0] += coef[is] * scale[2] * src[m_idx(0, 0, is)];
                }
            }
        }
    }
    //------ X +
    for (lid_t i2 = 0; i2 < M_N; i2++) {
        for (lid_t i1 = 0; i1 < M_N; i1++) {
            for (lid_t i0 = M_N - M_GS; i0 < M_N; i0++) {
                // get the data pointer in front of the row for every cache line
                real_p       trg = data_trg + m_idx(i0, i1, i2);  // cache line for writting
                const real_p src = data_src + m_idx(i0, i1, i2);  // cache line for reading
                // loop on the stencil block
                trg[0] = 0.0;
                for (int is = -(length / 2); is <= (length / 2); is++) {
                    trg[0] += coef[is] * scale[0] * src[m_idx(is, 0, 0)];
                    trg[0] += coef[is] * scale[1] * src[m_idx(0, is, 0)];
                    trg[0] += coef[is] * scale[2] * src[m_idx(0, 0, is)];
                }
            }
        }
    }
    //------ Y -
    for (lid_t i2 = 0; i2 < M_N; i2++) {
        for (lid_t i1 = 0; i1 < M_GS; i1++) {
            for (lid_t i0 = 0; i0 < M_N; i0++) {
                // get the data pointer in front of the row for every cache line
                real_p       trg = data_trg + m_idx(i0, i1, i2);  // cache line for writting
                const real_p src = data_src + m_idx(i0, i1, i2);  // cache line for reading
                // loop on the stencil block
                trg[0] = 0.0;
                for (int is = -(length / 2); is <= (length / 2); is++) {
                    trg[0] += coef[is] * scale[0] * src[m_idx(is, 0, 0)];
                    trg[0] += coef[is] * scale[1] * src[m_idx(0, is, 0)];
                    trg[0] += coef[is] * scale[2] * src[m_idx(0, 0, is)];
                }
            }
        }
    }
    //------ Y +
    for (lid_t i2 = 0; i2 < M_N; i2++) {
        for (lid_t i1 = M_N - M_GS; i1 < M_N; i1++) {
            for (lid_t i0 = 0; i0 < M_N; i0++) {
                // get the data pointer in front of the row for every cache line
                real_p       trg = data_trg + m_idx(i0, i1, i2);  // cache line for writting
                const real_p src = data_src + m_idx(i0, i1, i2);  // cache line for reading
                // loop on the stencil block
                trg[0] = 0.0;
                for (int is = -(length / 2); is <= (length / 2); is++) {
                    trg[0] += coef[is] * scale[0] * src[m_idx(is, 0, 0)];
                    trg[0] += coef[is] * scale[1] * src[m_idx(0, is, 0)];
                    trg[0] += coef[is] * scale[2] * src[m_idx(0, 0, is)];
                }
            }
        }
    }
    //------ Z -
    for (lid_t i2 = 0; i2 < M_GS; i2++) {
        for (lid_t i1 = 0; i1 < M_N; i1++) {
            for (lid_t i0 = 0; i0 < M_N; i0++) {
                // get the data pointer in front of the row for every cache line
                real_p       trg = data_trg + m_idx(i0, i1, i2);  // cache line for writting
                const real_p src = data_src + m_idx(i0, i1, i2);  // cache line for reading
                // loop on the stencil block
                trg[0] = 0.0;
                for (int is = -(length / 2); is <= (length / 2); is++) {
                    trg[0] += coef[is] * scale[0] * src[m_idx(is, 0, 0)];
                    trg[0] += coef[is] * scale[1] * src[m_idx(0, is, 0)];
                    trg[0] += coef[is] * scale[2] * src[m_idx(0, 0, is)];
                }
            }
        }
    }
    //------ Z +
    for (lid_t i2 = M_N - M_GS; i2 < M_N; i2++) {
        for (lid_t i1 = 0; i1 < M_N; i1++) {
            for (lid_t i0 = 0; i0 < M_N; i0++) {
                // get the data pointer in front of the row for every cache line
                real_p       trg = data_trg + m_idx(i0, i1, i2);  // cache line for writting
                const real_p src = data_src + m_idx(i0, i1, i2);  // cache line for reading
                // loop on the stencil block
                trg[0] = 0.0;
                for (int is = -(length / 2); is <= (length / 2); is++) {
                    trg[0] += coef[is] * scale[0] * src[m_idx(is, 0, 0)];
                    trg[0] += coef[is] * scale[1] * src[m_idx(0, is, 0)];
                    trg[0] += coef[is] * scale[2] * src[m_idx(0, 0, is)];
                }
            }
        }
    }
    // -------------------------------------------------------------------------
}

#endif  // SRC_LAPLACIAN_IPP_
