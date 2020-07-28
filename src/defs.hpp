#ifndef SRC_DEFS_HPP_
#define SRC_DEFS_HPP_

#include <execinfo.h>
#include <mpi.h>
#include <omp.h>
#include <p8est.h>

#include <cstdio>
#include <cstdlib>

/**
 * @name user changeable parameters 
 * @{
 */
#define M_N 16          //!< size of one block (M_N x M_N x M_N)
#define M_GS 4          //!< number of ghost points
#define M_ALIGNMENT 16  //!< memory alignement (in Byte, 16 = 2 doubles = 4 floats)

#define M_WAVELET_N 4
#define M_WAVELET_NT 2
/** @} */

/**
 * @name memory sizes shortcuts
 * @{
 */
#define M_HN (M_N / 2)  //!< half size of a block
#define M_STRIDE (2 * M_GS + M_N)
/** @} */

/**
 * @name memory management
 * @{
 */
#define M_MPI_REAL MPI_DOUBLE  //!< type used for the MPI communication (double by default)

/**
 * @brief returns true if the memory is aligned
 * 
 */
#define m_isaligned(a)                      \
    ({                                      \
        const void* a_ = (void*)(a);        \
        ((uintptr_t)a_) % M_ALIGNMENT == 0; \
    })

#if defined(__INTEL_COMPILER)
#define m_assume_aligned(a)                                  \
    ({                                                       \
        __typeof__(a) a_ = (a);                              \
        m_assert(m_isaligned(a_), "data has to be aligned"); \
        __assume_aligned(a_, M_ALIGNMENT);                   \
    })
#define m_calloc(size)                                 \
    ({                                                 \
        size_t size_ = (size_t)(size);                 \
        void*  data  = _mm_malloc(size_, M_ALIGNMENT); \
        memset(data, 0, size_);                        \
        data;                                          \
    })
#define m_free(data)                     \
    ({                                   \
        __typeof__(data) data_ = (data); \
        _mm_free((void*)data_);          \
    })
#else  //defined(__GNUC__)
#define m_assume_aligned(a)                                  \
    ({                                                       \
        __typeof__(a) a_ = (a);                              \
        m_assert(m_isaligned(a_), "data has to be aligned"); \
        __builtin_assume_aligned(a_, M_ALIGNMENT);           \
    })
/**
 * @brief allocate a given size (in Byte) and set to 0 the array.
 * the return pointer is aligned to M_ALIGMEMENT
 */
#define m_calloc(size)                                    \
    ({                                                    \
        size_t size_ = (size_t)(size);                    \
        void*  data  = aligned_alloc(M_ALIGNMENT, size_); \
        memset(data, 0, size_);                           \
        data;                                             \
    })
/**
 * @brief frees the pointer allocated using @ref m_calloc()
 */
#define m_free(data)                     \
    ({                                   \
        __typeof__(data) data_ = (data); \
        free((void*)data_);              \
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
#define m_max(a, b)             \
    ({                          \
        __typeof__(a) a_ = (a); \
        __typeof__(b) b_ = (b); \
        a_ > b_ ? a_ : b_;      \
    })

/**
 * @brief returns the min of two expressions
 * 
 */
#define m_min(a, b)             \
    ({                          \
        __typeof__(a) a_ = (a); \
        __typeof__(b) b_ = (b); \
        a_ < b_ ? a_ : b_;      \
    })

/**
 * @brief returns the sign of a number
 * 
 */
#define m_sign(a)                    \
    ({                               \
        __typeof__(a) a_    = (a);   \
        __typeof__(a) zero_ = 0;     \
        (zero_ < a_) - (a_ < zero_); \
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
#define m_pos(pos, i0, i1, i2, hgrid, xyz) \
    ({                                     \
        __typeof__(i0) i0_ = (i0);         \
        __typeof__(i1) i1_ = (i1);         \
        __typeof__(i2) i2_ = (i2);         \
                                           \
        pos[0] = i0_ * hgrid[0] + xyz[0];  \
        pos[1] = i1_ * hgrid[1] + xyz[1];  \
        pos[2] = i2_ * hgrid[2] + xyz[2];  \
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
#define m_pos_relative(offset, i0, i1, i2, hgrid) \
    ({                                            \
        __typeof__(i0) i0_ = (i0);                \
        __typeof__(i1) i1_ = (i1);                \
        __typeof__(i2) i2_ = (i2);                \
                                                  \
        offset[0] = i0_ * hgrid[0];               \
        offset[1] = i1_ * hgrid[1];               \
        offset[2] = i2_ * hgrid[2];               \
    })

/**
 * @brief returns the size (in elements) of one block
 * 
 */
#define m_blockmemsize(lda)                              \
    ({                                                   \
        __typeof__(lda) lda_ = (lda);                    \
        (size_t)(lda_ * M_STRIDE * M_STRIDE * M_STRIDE); \
    })

/**
 * @brief returns the memory offset of the first block element: (0,0,0)
 * 
 */
#define m_zeroidx(ida, mem)                                        \
    ({                                                             \
        __typeof__(mem) mem_ = (mem);                              \
        sid_t  ida_          = (ida);                              \
        lid_t  gs_           = (mem_->gs());                       \
        size_t str_          = (size_t)(mem_->stride());           \
        (size_t)(gs_ + str_ * (gs_ + str_ * (gs_ + str_ * ida_))); \
    })

/**
 * @brief returns the memory offset to reach a 3D position (i0,i1,i2) given a stride (str).
 * 
 * This macro is to be used with the function GridBlock::data()
 * 
 * @note: we cast the stride to size_t to ensure a proper conversion while computing the adress
 */
#define m_sidx(i0, i1, i2, ida, str)                               \
    ({                                                             \
        lid_t  i0_  = (i0);                                        \
        lid_t  i1_  = (i1);                                        \
        lid_t  i2_  = (i2);                                        \
        sid_t  ida_ = (ida);                                       \
        size_t str_ = (size_t)(str);                               \
        (size_t)(i0_ + str_ * (i1_ + str_ * (i2_ + str_ * ida_))); \
    })

/**
 * @brief returns the memory offset to reach a 3D position (i0,i1,i2) given a @ref MemLayout
 * 
 * This macro is equivalent to @ref m_sidx with a stride given by MemLayout::stride()
 * 
 */
#define m_midx(i0, i1, i2, ida, mem)            \
    ({                                          \
        m_sidx(i0, i1, i2, ida, mem->stride()); \
    })

/**
 * @brief returns the memory index given 3D position, for a GridBlock object (only!)
 * 
 */
#define m_idx(i0, i1, i2)                \
    ({                                   \
        m_sidx(i0, i1, i2, 0, M_STRIDE); \
    })

/**
 * @brief returns the lenght of a quadrant at a given level,
 * assuming one octree is a cubic domain: (1 x 1 x 1)
 */
#define m_quad_len(level)                                  \
    ({                                                     \
        __typeof__(level) lvl_ = (level);                  \
        1.0 / (P8EST_ROOT_LEN / P8EST_QUADRANT_LEN(lvl_)); \
    })

/** @} */

/**
 * @name logs and verbosity 
 * 
 */
/**
 * @brief m_log will be displayed as a log, either by every rank or only by the master (given LOG_ALLRANKS)
 * 
 */
#ifndef LOG_MUTE
#ifndef LOG_ALLRANKS
#define m_log(format, ...)                                 \
    ({                                                     \
        int rank;                                          \
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);              \
        if (rank == 0) {                                   \
            char def_nhyipns[1024];                        \
            sprintf(def_nhyipns, format, ##__VA_ARGS__);   \
            fprintf(stdout, "[murphy] %s\n", def_nhyipns); \
        }                                                  \
    })
#else
#define m_log(format, ...)                                      \
    ({                                                          \
        int rank;                                               \
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);                   \
        char def_nhyipns[1024];                                 \
        sprintf(def_nhyipns, format, ##__VA_ARGS__);            \
        fprintf(stdout, "[%d murphy] %s\n", rank, def_nhyipns); \
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
#define m_verb(format, ...)                                     \
    ({                                                          \
        int rank;                                               \
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);                   \
        char def_nhyipns[1024];                                 \
        sprintf(def_nhyipns, format, ##__VA_ARGS__);            \
        fprintf(stdout, "[%d murphy] %s\n", rank, def_nhyipns); \
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
#define m_assert(cond, ...)                                                                                            \
    ({                                                                                                                 \
        if (!(cond)) {                                                                                                 \
            char def_nhyzpns[1024];                                                                                    \
            sprintf(def_nhyzpns, __VA_ARGS__);                                                                         \
            fprintf(stdout, "[murphy-assert] '%s' FAILED: %s at %s (l:%d)\n", #cond, def_nhyzpns, __FILE__, __LINE__); \
            fflush(stdout);                                                                                            \
            MPI_Abort(MPI_COMM_WORLD, MPI_ERR_ASSERT);                                                                 \
        }                                                                                                              \
    })
#endif

/**
 * @brief entry and exit of functions, enabled if VERBOSE is enabled
 * 
 */
#ifdef VERBOSE
#define m_begin                                                                             \
    m_assert(omp_get_num_threads() == 1, "no MPI is allowed in an openmp parallel region"); \
    double m_def_T0 = MPI_Wtime();                                                          \
    m_verb("----- entering %s", __func__);
#define m_end                                                                               \
    m_assert(omp_get_num_threads() == 1, "no MPI is allowed in an openmp parallel region"); \
    double m_def_T1 = MPI_Wtime();                                                          \
    m_verb("----- leaving %s after %lf [s]", __func__, (m_def_T1) - (m_def_T0));
#else
#define m_begin \
    { ((void)0); }
#define m_end \
    { ((void)0); }
#endif
/** @} */

#endif  // SRC_DEFS_HPP_
