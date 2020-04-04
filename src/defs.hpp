#ifndef SRC_DEFS_HPP_
#define SRC_DEFS_HPP_

#include <mpi.h>
#include <cstdio>

#define M_N 16
#define M_GS 2
#define M_ALIGNMENT 16  // the alignement in bytes
#define M_MPI_REAL MPI_DOUBLE

#define M_DN 2 * MN  // the double size of one block
#define M_HN MN / 2  // the half size of one block
#define M_STRIDE (2 * M_GS + M_N)

#if defined(__INTEL_COMPILER)
#define m_assume_aligned(a)          \
    ({                               \
        __typeof__(a) a_ = (a);      \
        __assume_aligned(a_, FLUPS); \
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
#define m_assume_aligned(a)                  \
    ({                                       \
        __typeof__(a) a_ = (a);              \
        __builtin_assume_aligned(a_, FLUPS); \
    })
#define m_calloc(size)                                            \
    ({                                                            \
        void*  data;                                              \
        size_t size_ = (size_t)(size);                            \
        int    err   = posix_memalign(&data, M_ALIGNMENT, size_); \
        memset(data, 0, size_);                                   \
        data;                                                     \
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

#define m_blockmemsize(lda)                                \
    ({                                                  \
        __typeof__(lda) _lda = (lda);                   \
        (size_t)(lda * M_STRIDE * M_STRIDE * M_STRIDE); \
    })

/**
 * @brief return a memory index given 3D position (i0,i1,i2), a stride (str) and a ghost size (gs)
 * 
 */
#define m_sidx(i0, i1, i2, ida, str, gs)                                       \
    ({                                                                         \
        __typeof__(i0) _i0   = (i0);                                           \
        __typeof__(i1) _i1   = (i1);                                           \
        __typeof__(i2) _i2   = (i2);                                           \
        __typeof__(ida) _ida = (ida);                                          \
        __typeof__(str) _str = (str);                                          \
        __typeof__(gs) _gs   = (gs);                                           \
        (_gs + _i0) + _str*((_gs + _i1) + _str * ((_gs + _i2) + _str * _ida)); \
    })

/**
 * @brief return the memory index given 3D position and a dimension number
 * 
 */
#define m_idx(i0, i1, i2, ida)                   \
    ({                                           \
        m_sidx(i0, i1, i2, ida, M_STRIDE, M_GS); \
    })

/**
 * @brief m_log will be displayed as a log
 * 
 */
#ifndef LOG_ALLRANKS
#define m_log(format, ...)                               \
    ({                                                   \
        int rank;                                        \
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);            \
        if (rank == 0) {                                 \
            char def_nhyipns[1024];                      \
            sprintf(def_nhyipns, format, ##__VA_ARGS__); \
            fprintf(stdout, "[murphy] %s\n", def_nhyipns); \
        }                                                \
    })
#else
#define m_log(format, ...)                                    \
    ({                                                        \
        int rank;                                             \
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);                 \
        char def_nhyipns[1024];                               \
        sprintf(def_nhyipns, format, ##__VA_ARGS__);          \
        fprintf(stdout, "[%d murphy] %s\n", rank, def_nhyipns); \
    })
#endif
/**
 * @brief m_verb will be displayed if VERBOSE is enabled
 * 
 */
#ifdef VERBOSE
#define m_verb(format, ...)                                   \
    ({                                                        \
        int rank;                                             \
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);                 \
        char def_nhyipns[1024];                               \
        sprintf(def_nhyipns, format, ##__VA_ARGS__);          \
        fprintf(stdout, "[%d murphy] %s\n", rank, def_nhyipns); \
    })
#else
#define m_verb(format, ...) \
    { ((void)0); }
#endif

/**
 * @brief entry and exit of functions, enabled if VERBOSE is enabled
 * 
 */
#define m_begin                         \
    double def_idajfl_T0 = MPI_Wtime(); \
    m_verb("----- entering %s", __func__);

#define m_end                           \
    double def_idajfl_T1 = MPI_Wtime(); \
    m_verb("----- leaving %s after %lf [s]", __func__, (def_idajfl_T1) - (def_idajfl_T0));

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
            fprintf(stderr, "[murphy-assert] '%s' FAILED: %s at %s (l:%d)", #cond, def_nhyzpns, __FILE__, __LINE__); \
            MPI_Abort(MPI_COMM_WORLD, MPI_ERR_ASSERT);                                                               \
        }                                                                                                            \
    })
#endif

#endif  // SRC_DEFS_HPP_
