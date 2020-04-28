#ifndef SRC_MURPHY_HPP_
#define SRC_MURPHY_HPP_

#include <p8est.h>

#include <map>
#include <string>

#include "defs.hpp"

typedef int32_t lid_t;
typedef int8_t  sid_t;
typedef double  real_t;

typedef real_t* __restrict real_p;

typedef std::map<std::string, sid_t>  ldamap_t;
typedef std::map<std::string, real_p> datamap_t;

typedef p8est_quadrant_t qdrt_t;

/**
 * @brief id of a quadrant
 * 
 */
typedef struct qid_t {
    p4est_locidx_t cid;  //!< cummulative_id, does not depends on the tree id
    p4est_locidx_t qid;  //!< quadrant id = depends on the tree
    p4est_topidx_t tid;  //!< tree id
} qid_t;

/**
 * @brief defines the different supported boundary conditions
 * 
 */
typedef enum bctype_t {
    M_BC_EVEN,   //!< EVEN condition: an EVEN bondary condition with respect to a given flux at the interface
    M_BC_ODD,    //!< ODD condition: imposes a zero value wrt the interface
    M_BC_ZERO,   //!< set 0 outside the computational domain
    M_BC_EXTRAP  //!< extrapolate outside the computational domain
} bctype_t;

#endif  // SRC_MURPHY_HPP_
