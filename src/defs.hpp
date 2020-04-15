#ifndef SRC_DEFS_HPP_
#define SRC_DEFS_HPP_

#include <mpi.h>
#include <omp.h>
#include <cstdio>
#include <p8est.h>

#include <cstdlib>

#define M_N 16
#define M_GS 2
#define M_ALIGNMENT 16  // the alignement in bytes
#define M_MPI_REAL MPI_DOUBLE

#define M_DN (2 * M_N)  // the double size of one block
#define M_HN (M_N / 2)  // the half size of one block
#define M_STRIDE (2 * M_GS + M_N)

#define m_isaligned(a)                      \
    ({                                      \
        const void* _a = (a);               \
        ((uintptr_t)_a) % M_ALIGNMENT == 0; \
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
#define m_free(data)    \
    ({                  \
        _mm_free(data); \
    })
#elif defined(__GNUC__)
#define m_assume_aligned(a)                                  \
    ({                                                       \
        __typeof__(a) a_ = (a);                              \
        m_assert(m_isaligned(a_), "data has to be aligned"); \
        __builtin_assume_aligned(a_, M_ALIGNMENT);           \
    })
#define m_calloc(size)                                    \
    ({                                                    \
        size_t size_ = (size_t)(size);                    \
        void*  data  = aligned_alloc(M_ALIGNMENT, size_); \
        memset(data, 0, size_);                           \
        data;                                             \
    })
#define m_free(data) \
    ({               \
        free(data);  \
    })
#endif

#define m_max(a, b)             \
    ({                          \
        __typeof__(a) _a = (a); \
        __typeof__(b) _b = (b); \
        _a > _b ? _a : _b;      \
    })

#define m_min(a, b)             \
    ({                          \
        __typeof__(a) _a = (a); \
        __typeof__(b) _b = (b); \
        _a < _b ? _a : _b;      \
    })

#define m_pos(pos, i0, i1, i2, hgrid, xyz)        \
    ({                                            \
        __typeof__(i0) _i0 = (i0);                \
        __typeof__(i1) _i1 = (i1);                \
        __typeof__(i2) _i2 = (i2);                \
                                                  \
        pos[0] = (_i0 + 0.5) * hgrid[0] + xyz[0]; \
        pos[1] = (_i1 + 0.5) * hgrid[1] + xyz[1]; \
        pos[2] = (_i2 + 0.5) * hgrid[2] + xyz[2]; \
    })

#define m_blockmemsize(lda)                             \
    ({                                                  \
        __typeof__(lda) _lda = (lda);                   \
        (size_t)(lda * M_STRIDE * M_STRIDE * M_STRIDE); \
    })

/**
 * @brief returns the memory position of (0,0,0)
 * 
 */
#define m_zeroidx(ida, mem)                                        \
    ({                                                             \
        sid_t  ida_ = (ida);                                       \
        lid_t  gs_  = (mem->gs());                                 \
        size_t str_ = (size_t)(mem->stride());                     \
        (size_t)(gs_ + str_ * (gs_ + str_ * (gs_ + str_ * ida_))); \
    })

/**
 * @brief return the shift in memory to reach a 3D position (i0,i1,i2) given a stride (str).
 * 
 * This macro is to be used with the function GridBlock::data()
 * 
 * @note: we cast the stride to size_t to ensure a proper conversion while computing the adress
 * 
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
 * @brief return the shift in memory to reach a 3D position (i0,i1,i2) given a MemLayout mem
 * 
 * This macro is equivalent to @ref m_sidx with a stride given by MemLayout::stride()
 * 
 */
#define m_midx(i0, i1, i2, ida, mem)            \
    ({                                          \
        m_sidx(i0, i1, i2, ida, mem->stride()); \
    })

/**
 * @brief return the memory index given 3D position and a dimension number
 * 
 */
#define m_idx(i0, i1, i2)                \
    ({                                   \
        m_sidx(i0, i1, i2, 0, M_STRIDE); \
    })

/**
 * @brief return the lenght of a quadrant at a given level
 * 
 */
#define m_quad_len(level)                                  \
    ({                                                     \
        __typeof__(level) lvl_ = (level);                  \
        1.0 / (P8EST_ROOT_LEN / P8EST_QUADRANT_LEN(lvl_)); \
    })

/**
 * @brief m_log will be displayed as a log
 * 
 */
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
#define m_assert(cond, ...)                                                                                          \
    ({                                                                                                               \
        if (!(cond)) {                                                                                               \
            char def_nhyzpns[1024];                                                                                  \
            sprintf(def_nhyzpns, __VA_ARGS__);                                                                       \
            fprintf(stderr, "[murphy-assert] '%s' FAILED: %s at %s (l:%d)\n", #cond, def_nhyzpns, __FILE__, __LINE__); \
            MPI_Abort(MPI_COMM_WORLD, MPI_ERR_ASSERT);                                                               \
        }                                                                                                            \
    })
#endif


/**
 * @brief entry and exit of functions, enabled if VERBOSE is enabled
 * 
 */
#define m_begin                                                                            \
    m_assert(omp_get_num_threads() == 1, "no MPI is allowed in an openmp parallel region"); \
    double def_idajfl_T0 = MPI_Wtime();                                                    \
    m_verb("----- entering %s", __func__);

#define m_end                                                                              \
    m_assert(omp_get_num_threads() == 1, "no MPI is allowed in an openmp parallel region"); \
    double def_idajfl_T1 = MPI_Wtime();                                                    \
    m_verb("----- leaving %s after %lf [s]", __func__, (def_idajfl_T1) - (def_idajfl_T0));

#endif  // SRC_DEFS_HPP_