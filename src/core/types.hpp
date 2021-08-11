#ifndef SRC_CORE_TYPES_HPP_
#define SRC_CORE_TYPES_HPP_

#include <p8est.h>
#include <functional>

//------------------------------------------------------------------------------
/**
 * @name typename for standard variables
 * @{
 */
using iface_t  = short;             //!< face index type (0->26)
using lda_t    = short;             //!< leading dimension array type
using level_t  = short;             //!< data type for a level
using bidx_t   = int;               //!< indexing inside a block
using iter_t   = int;               //!< indexing iterations
using rank_t   = int;               //!< rank data type
using real_t   = double;            //!< floating point data
using iblock_t = p4est_locidx_t;    //!< datatype to index blocks, ghosts and quadrants
using qdrt_t   = p8est_quadrant_t;  //!< p4est_quadrant, the default type is too long
using short_t  = short;             //!< short integer, only used when nothing else fits: order of a method, etc...
using long_t   = long;              //!< long integer (for global quantities)
/**@}*/

//-------------------------------------------------------------------------
template<typename R,typename... T>
using lambda_t = std::function<R(T...)>;

/**
 * @brief definition of the lambda functions
 */
template <typename R, typename... T>
using lambda_i3_t = lambda_t<R, const bidx_t, const bidx_t, const bidx_t, T...>;
// using lambda_i3block_t = std::function<R(const bidx_t i0, const bidx_t i1, const bidx_t i2, T...)>;

// defines a function that takes 3 arguments and returns an offset
using accessor_t = std::function<bidx_t(const bidx_t, const bidx_t, const bidx_t, const bidx_t)>;

//------------------------------------------------------------------------------
/**
 * @brief id of a quadrant
 * 
 */
typedef struct qid_t {
    p4est_locidx_t cid;  //!< cummulative_id = the absolute id of the quadrant, used to access arrays
    p4est_topidx_t tid;  //!< tree id
    // actual id of the quadrant: depend if it's a quadrant, a mirror or a ghost
    union {
        p4est_locidx_t qid;  //!< if a normal quadrant: quadrant id = depends on the tree
        p4est_locidx_t mid;  //!< if a mirror quadrant: mirror id
        // p4est_locidx_t gid;  //!< if a ghost quadrant: ghost id
    };
} qid_t;

//------------------------------------------------------------------------------
/**
 * @brief defines the different supported boundary conditions, the order is driven by M_WAVELET_N
 */
typedef enum bctype_t {
    M_BC_NONE,    //!< no boundary condition is given
    M_BC_NEU,     //!< EVEN condition: an EVEN bondary condition with respect to a given flux at the interface
    M_BC_DIR,     //!< ODD condition: imposes a zero value wrt the interface
    M_BC_ZERO,    //!< set 0 outside the computational domain
    M_BC_EXTRAP,  //!< extrapolate outside the computational domain with polynomial
} bctype_t;

typedef enum m_direction_t {
    M_FORWARD,
    M_BACKWARD
} m_direction_t;

//-------------------------------------------------------------------------
// must be removed from the code
typedef int  lid_t;
typedef char sid_t;
typedef real_t* __restrict real_p;  //!< pointer type = root of the memory allocation

#endif  // SRC_CORE_TYPES_HPP_
