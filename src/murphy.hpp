#ifndef SRC_MURPHY_HPP_
#define SRC_MURPHY_HPP_

#include <p8est.h>

#include <map>
#include <string>
#include <limits>

#include "defs.hpp"

/**
 * @name typename for indexes
 * @{
 */
using iface_t  = int8_t;          //!< face index type (0->26)
using lda_t    = int8_t;          //!< leading dimension array type
using level_t  = int8_t;          //!< data type for a level
using rank_t   = int;             //!< rank data type
using iblock_t = p4est_locidx_t;  //!< datatype to index blocks, ghosts and quadrants
using qdrt_t   = p8est_quadrant_t;
/**@}*/

/**
 * @name typename for data pointers
 * @{
 */
using real_t        = double;                    //!< data subtype
using data_ptr      = real_t* __restrict__;        //!< data pointer type = root of the data, i.e. point (0,0,0)
using mem_ptr       = real_t* __restrict__;        //!< pointer type = root of the memory allocation
using const_mem_ptr = const real_t* __restrict__;  //!< pointer type = root of the memory allocation with the const
/**@}*/

// must be removed from the code

typedef int32_t lid_t;
typedef int8_t  sid_t;
typedef real_t* __restrict real_p;  //!< pointer type = root of the memory allocation
typedef std::map<std::string, real_p> datamap_t;

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

/**
 * @brief defines the different supported boundary conditions, the order is driven by M_WAVELET_N
 */
typedef enum bctype_t {
    M_BC_NONE,    //!< no boundary condition is given
    M_BC_NEU,    //!< EVEN condition: an EVEN bondary condition with respect to a given flux at the interface
    M_BC_DIR,     //!< ODD condition: imposes a zero value wrt the interface
    M_BC_ZERO,    //!< set 0 outside the computational domain
    M_BC_EXTRAP,  //!< extrapolate outside the computational domain with polynomial
} bctype_t;

typedef enum m_direction_t {
    M_FORWARD,
    M_BACKWARD
} m_direction_t;

// template <typename T>
// class m_shared_ptr {
//    public:
//     T* ptr;
//     T* operator()() { return ptr; };
// }

// template <typename T>
// class m_owned_ptr:public m_shared_ptr{};

void murphy_init(int argc, char* argv[]);
void murphy_finalize();

#endif  // SRC_MURPHY_HPP_
