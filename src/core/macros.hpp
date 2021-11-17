#ifndef SRC_CORE_MACRO_HPP_
#define SRC_CORE_MACRO_HPP_

#include <execinfo.h>
#include <mpi.h>
#include <omp.h>
#include <p8est.h>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>

#include "core/types.hpp"
#include "hdf5.h"

#ifndef MPI_NONASYNC
#define M_MPI_AGGRESSIVE
#endif

// register the current git commit for tracking purpose
#ifdef GIT_COMMIT
#define M_GIT_COMMIT GIT_COMMIT
#else
#define M_GIT_COMMIT "?"
#endif

// check if the compilation defines the order of the wavelet. if not, we do it
#ifndef WAVELET_N
#define M_WAVELET_N 2
#else
#define M_WAVELET_N WAVELET_N
#endif

#ifndef WAVELET_NT
#define M_WAVELET_NT 0
#else
#define M_WAVELET_NT WAVELET_NT
#endif

/**
 * @name user-defined parameters 
 * @{
 */
#define M_N 24          //!< size of one block (M_N x M_N x M_N), must be EVEN
#define M_ALIGNMENT 16  //!< memory alignement (in Byte, 16 = 2 doubles = 4 floats)

#ifndef BLOCK_GS
#define M_GS 4  //!< memory space for the ghost points (not the actual number of ghost points!)
#else
#define M_GS BLOCK_GS
#endif
/** @} */

// /**
//  * @name memory sizes shortcuts
//  * @{
//  */
#define M_NCENTER (M_N / 2)        //!< center ID of a block
#define M_NHALF (M_N / 2)          //!< number of points in a half block
#define M_STRIDE (2 * M_GS + M_N)  //!< stride in memory

// // #ifndef BLOCK_GS
// // static constexpr bidx_t m_gs = 4;
// // #else
// // static constexpr bidx_t m_gs = BLOCK_GS;
// // #endif

// // static constexpr bidx_t m_stride_reg = (2 * m_gs + M_N);
// // static constexpr bidx_t m_stride_ext = ((m_stride_reg * sizeof(real_t)) % M_ALIGNMENT == 0) ?
// // static constexpr bidx_t m_stride[3] = {() % (M_ALIGNMENT / sizeof(real_t))}
// /** @} */


#define M_INLINE __attribute__((always_inline)) inline

/**
 * @name memory management
 * @{
 */
#define M_MPI_REAL MPI_DOUBLE          //!< type used for the MPI communication (double by default)
#define M_HDF5_REAL H5T_NATIVE_DOUBLE  //!< type used for the MPI communication (double by default)

/**
 * @brief returns true if the memory is aligned
 * 
 */
#define m_isaligned(a)                                  \
    ({                                                  \
        const void* m_isaligned_a_ = (void*)(a);        \
        ((uintptr_t)m_isaligned_a_) % M_ALIGNMENT == 0; \
    })

#if defined(__INTEL_COMPILER)
#define m_assume_aligned(a)                                                   \
    ({                                                                        \
        __typeof__(a) m_assume_aligned_a_ = (a);                              \
        m_assert(m_isaligned(m_assume_aligned_a_), "data has to be aligned"); \
        __assume_aligned(m_assume_aligned_a_, M_ALIGNMENT);                   \
    })
#define m_calloc(size)                                                                    \
    ({                                                                                    \
        size_t m_calloc_size_        = (size_t)(size) + M_ALIGNMENT - 1;                  \
        size_t m_calloc_padded_size_ = (m_calloc_size_) - (m_calloc_size_ % M_ALIGNMENT); \
        void*  m_calloc_data_        = _mm_malloc(m_calloc_padded_size_, M_ALIGNMENT);    \
        std::memset(m_calloc_data_, 0, m_calloc_padded_size_);                            \
        m_calloc_data_;                                                                   \
    })
#define m_free(data)                        \
    ({                                      \
        void* m_free_data_ = (void*)(data); \
        _mm_free(m_free_data_);             \
    })
#else  //defined(__GNUC__)
#define m_assume_aligned(a)                                                   \
    ({                                                                        \
        __typeof__(a) m_assume_aligned_a_ = (a);                              \
        m_assert(m_isaligned(m_assume_aligned_a_), "data has to be aligned"); \
        __builtin_assume_aligned(m_assume_aligned_a_, M_ALIGNMENT);           \
    })
/**
 * @brief allocate a given size (in Byte) and set to 0 the array.
 * the return pointer is aligned to M_ALIGMEMENT
 */
#define m_calloc(size)                                                                    \
    ({                                                                                    \
        size_t m_calloc_size_        = (size_t)(size) + M_ALIGNMENT - 1;                  \
        size_t m_calloc_padded_size_ = (m_calloc_size_) - (m_calloc_size_ % M_ALIGNMENT); \
        void * m_calloc_data_;\
        posix_memalign(&m_calloc_data_,M_ALIGNMENT,m_calloc_padded_size_);\
        std::memset(m_calloc_data_, 0, m_calloc_padded_size_);                            \
        m_calloc_data_;                                                                   \
    })
/**
 * @brief frees the pointer allocated using @ref m_calloc()
 */
#define m_free(data)                        \
    ({                                      \
        void* m_free_data_ = (void*)(data); \
        std::free(m_free_data_);            \
    })
#endif
/** @} */

/**
 * @name min max sign macros
 * 
 */
/**
 * @brief returns the max of two expressions
 * 
 */
#define m_max(a, b)                                      \
    ({                                                   \
        __typeof__(a) m_max_a_ = (a);                    \
        __typeof__(b) m_max_b_ = (b);                    \
        (m_max_a_ > m_max_b_) ? (m_max_a_) : (m_max_b_); \
    })

/**
 * @brief returns the min of two expressions
 * 
 */
#define m_min(a, b)                                      \
    ({                                                   \
        __typeof__(a) m_min_a_ = (a);                    \
        __typeof__(b) m_min_b_ = (b);                    \
        (m_min_a_ < m_min_b_) ? (m_min_a_) : (m_min_b_); \
    })

/**
 * @brief returns the sign of a number: +1 if positive, 0 if 0 and -1 if negative
 * 
 */
#define m_sign(a)                                                \
    ({                                                           \
        __typeof__(a) m_sign_a_    = (a);                        \
        __typeof__(a) m_sign_zero_ = 0;                          \
        (m_sign_zero_ < m_sign_a_) - (m_sign_a_ < m_sign_zero_); \
    })

#define m_fequal(a, b)                                                                 \
    ({                                                                                 \
        real_t m_equal_a_ = (a);                                                       \
        real_t m_equal_b_ = (b);                                                       \
        (std::fabs(m_equal_a_ - m_equal_b_) < std::numeric_limits<real_t>::epsilon()); \
    })

/** @} */

/**
 * @name Block related operations
 * 
 */
/**
 * @brief returns the position of a point (i0,i1,i2) wrt to the computational domain
 * 
 * @param i0 the index in the x direction within a block
 * @param i1 the index in the y direction within a block
 * @param i2 the index in the z direction within a block
 * @param hgrid the local grid spacing
 * @param xyz the position of the origin of the block (!= the position of (0,0,0))
 * 
 */
#define m_pos(pos, i0, i1, i2, hgrid, xyz)      \
    ({                                          \
        __typeof__(i0) m_pos_i0_ = (i0);        \
        __typeof__(i1) m_pos_i1_ = (i1);        \
        __typeof__(i2) m_pos_i2_ = (i2);        \
                                                \
        pos[0] = m_pos_i0_ * hgrid[0] + xyz[0]; \
        pos[1] = m_pos_i1_ * hgrid[1] + xyz[1]; \
        pos[2] = m_pos_i2_ * hgrid[2] + xyz[2]; \
    })

/**
 * @brief returns the position of a point (i0,i1,i2) wrt ot the origin of the block
 * 
 * @param i0 the index in the x direction within a block
 * @param i1 the index in the y direction within a block
 * @param i2 the index in the z direction within a block
 * @param hgrid the local grid spacing
 * @param xyz the position of the left corner of the block
 * 
 */
#define m_pos_relative(offset, i0, i1, i2, hgrid)  \
    ({                                             \
        __typeof__(i0) m_pos_relative_i0_ = (i0);  \
        __typeof__(i1) m_pos_relative_i1_ = (i1);  \
        __typeof__(i2) m_pos_relative_i2_ = (i2);  \
                                                   \
        offset[0] = m_pos_relative_i0_ * hgrid[0]; \
        offset[1] = m_pos_relative_i1_ * hgrid[1]; \
        offset[2] = m_pos_relative_i2_ * hgrid[2]; \
    })

// /**
//  * @brief returns the size (in elements) of one block
//  *
//  */
// #define m_blockmemsize(lda)                                             \
//     ({                                                                  \
//         __typeof__(lda) m_blockmemsize_lda_ = (lda);                    \
//         (bidx_t)(m_blockmemsize_lda_ * M_STRIDE * M_STRIDE * M_STRIDE); \
//     })

/**
 * @brief returns the memory offset of the first block element: (0,0,0)
 * 
 */
#define m_zeroidx(ida, mem)                                                                                                              \
    ({                                                                                                                                   \
        __typeof__(mem) m_zeroidx_mem_ = (mem);                                                                                          \
        lda_t           m_zeroidx_ida_ = (lda_t)(ida);                                                                                   \
        bidx_t          m_zeroidx_gs_  = (bidx_t)(m_zeroidx_mem_->gs());                                                                 \
        bidx_t          m_zeroidx_str_ = (bidx_t)(m_zeroidx_mem_->stride());                                                             \
        (bidx_t)(m_zeroidx_gs_ + m_zeroidx_str_ * (m_zeroidx_gs_ + m_zeroidx_str_ * (m_zeroidx_gs_ + m_zeroidx_str_ * m_zeroidx_ida_))); \
    })

/**
 * @brief returns the memory offset to reach a 3D position (i0,i1,i2) given a stride (str).
 * 
 * The ghost point position is not taken into account here as we already have the (0,0,0) position with GridBlock::data().
 * 
 * @note: we cast the stride to size_t to ensure a proper conversion while computing the adress
 */
// #define m_sidx(i0, i1, i2, ida, str)                                                                                \
//     ({                                                                                                              \
//         lid_t  m_sidx_i0_  = (lid_t)(i0);                                                                           \
//         lid_t  m_sidx_i1_  = (lid_t)(i1);                                                                           \
//         lid_t  m_sidx_i2_  = (lid_t)(i2);                                                                           \
//         lda_t  m_sidx_ida_ = (lda_t)(ida);                                                                          \
//         size_t m_sidx_str_ = (size_t)(str);                                                                         \
//         (size_t)(m_sidx_i0_ + m_sidx_str_ * (m_sidx_i1_ + m_sidx_str_ * (m_sidx_i2_ + m_sidx_str_ * m_sidx_ida_))); \
//     })

/**
 * @brief returns the memory offset to reach a 3D position (i0,i1,i2) given a @ref MemLayout
 * 
 * This macro is equivalent to @ref m_sidx with a stride given by MemLayout::stride()
 * The ghost point position is not taken into account here as we already have the (0,0,0) position with GridBlock::data().
 * 
 */
// #define m_midx(i0, i1, i2, ida, mem)                    \
//     ({                                                  \
//         __typeof__(mem) m_midx_mem_ = (mem);            \
//         m_sidx(i0, i1, i2, ida, m_midx_mem_->stride()); \
//     })

/**
 * @brief returns the memory offset to reach a 3D position (i0,i1,i2) given a @ref MemLayout
 * 
 * The ghost point position is not taken into account here as we already have the (0,0,0) position with GridBlock::data().
 * 
 */
// #define m_idx(i0, i1, i2)                \
//     ({                                   \
//         m_sidx(i0, i1, i2, 0, M_STRIDE); \
//     })

/** @} */

/**
 * @name logs and verbosity 
 * 
 */
extern short m_log_level_counter;
extern char  m_log_level_prefix[32];

#define m_log_level_plus                                  \
    ({                                                    \
        m_log_level_counter += (m_log_level_counter < 5); \
                                                          \
        m_log_level_prefix[0] = '\0';                     \
        for (short i = 0; i < m_log_level_counter; ++i) { \
            strcat(m_log_level_prefix, "  ");             \
        }                                                 \
    })
#define m_log_level_minus                                 \
    ({                                                    \
        m_log_level_counter -= (m_log_level_counter > 0); \
                                                          \
        m_log_level_prefix[0] = '\0';                     \
        for (short i = 0; i < m_log_level_counter; ++i) { \
            strcat(m_log_level_prefix, "  ");             \
        }                                                 \
    })
/**
 * @brief m_log will be displayed as a log, either by every rank or only by the master (given LOG_ALLRANKS)
 * 
 */
#ifndef LOG_MUTE
#ifndef LOG_ALLRANKS
#define m_log(format, ...)                                                       \
    ({                                                                           \
        int m_log_rank_;                                                         \
        MPI_Comm_rank(MPI_COMM_WORLD, &m_log_rank_);                             \
        if (m_log_rank_ == 0) {                                                  \
            char m_log_msg_[1024];                                               \
            sprintf(m_log_msg_, format, ##__VA_ARGS__);                          \
            fprintf(stdout, "[murphy] %s %s\n", m_log_level_prefix, m_log_msg_); \
        }                                                                        \
    })
#else
#define m_log(format, ...)                                                                   \
    ({                                                                                       \
        int m_log_rank_;                                                                     \
        MPI_Comm_rank(MPI_COMM_WORLD, &m_log_rank_);                                         \
        char m_log_msg_[1024];                                                               \
        sprintf(m_log_msg_, format, ##__VA_ARGS__);                                          \
        fprintf(stdout, "[%d murphy] %s %s\n", m_log_rank_, m_log_level_prefix, m_log_msg_); \
    })
#endif
#else
#define m_log(format, ...) \
    { ((void)0); }
#endif

/**
 * @brief m_verb will be displayed if VERBOSE is enabled
 * 
 */
#ifdef VERBOSE
#define m_verb(format, ...)                                             \
    ({                                                                  \
        int m_verb_rank_;                                               \
        MPI_Comm_rank(MPI_COMM_WORLD, &m_verb_rank_);                   \
        char m_verb_msg_[1024];                                         \
        sprintf(m_verb_msg_, format, ##__VA_ARGS__);                    \
        fprintf(stdout, "[%d murphy] %s\n", m_verb_rank_, m_verb_msg_); \
    })
#else
#define m_verb(format, ...) \
    { ((void)0); }
#endif

/**
 * @brief m_assert defines the assertion call, disable if NDEBUG is asked
 * 
 */
#ifdef NDEBUG
#define m_assert(cond, ...) \
    { ((void)0); }
#else
#define m_assert(cond, ...)                                                                                                               \
    ({                                                                                                                                    \
        bool m_assert_cond_ = (bool)(cond);                                                                                               \
        if (!(m_assert_cond_)) {                                                                                                          \
            char m_assert_msg_[1024];                                                                                                     \
            int  m_assert_rank_;                                                                                                          \
            MPI_Comm_rank(MPI_COMM_WORLD, &m_assert_rank_);                                                                               \
            sprintf(m_assert_msg_, __VA_ARGS__);                                                                                          \
            fprintf(stdout, "[%d murphy-assert] '%s' FAILED: %s (at %s:%d)\n", m_assert_rank_, #cond, m_assert_msg_, __FILE__, __LINE__); \
            fflush(stdout);                                                                                                               \
            MPI_Abort(MPI_COMM_WORLD, MPI_ERR_ASSERT);                                                                                    \
        }                                                                                                                                 \
    })
#endif

/**
 * @brief entry and exit of functions, enabled if VERBOSE is enabled
 * 
 */
#ifdef VERBOSE
#define m_begin                                                                             \
    m_assert(omp_get_num_threads() == 1, "no MPI is allowed in an openmp parallel region"); \
    double m_begin_T0 = MPI_Wtime();                                                        \
    m_verb("----- entering %s", __func__);
#define m_end                                                                               \
    m_assert(omp_get_num_threads() == 1, "no MPI is allowed in an openmp parallel region"); \
    double m_end_T1_ = MPI_Wtime();                                                         \
    m_verb("----- leaving %s after %lf [s]", __func__, (m_end_T1_) - (m_begin_T0));
#else
#define m_begin \
    { ((void)0); }
#define m_end \
    { ((void)0); }
#endif
/** @} */

#endif  // SRC_CORE_MACRO_HPP_
