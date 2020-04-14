#ifndef SRC_MURPHY_HPP_
#define SRC_MURPHY_HPP_

#include <p8est.h>
#include "defs.hpp"
#include <map>


typedef int32_t lid_t;
typedef int8_t  sid_t;
typedef double  real_t;

typedef real_t* __restrict real_p;

typedef std::map<std::string, sid_t> ldamap_t;
typedef std::map<std::string, real_p> datamap_t;

typedef p8est_quadrant_t qdrt_t;

// /**
//  * @brief definitions of the different errors
//  *
//  */
// typedef enum err_t {
//     M_SUCCESS,
//     M_ERROR,
//     M_FATAL
// } err_t;

/**
 * @brief id of a quadrant
 * 
 */
typedef struct qid_t {
    p4est_locidx_t cid; /*<cummulative_id, does not depends on the tree id */
    p4est_locidx_t qid; /*< quadrant id = depends on the tree */
    p4est_topidx_t tid; /*< tree id */
} qid_t;

/**
 * @brief extends the qid_t by adding a local ID
 * 
 * typicall used for ghosts
 * 
 */
typedef struct qid_ext_t : qid_t {
    p4est_locidx_t lid; /*<local id, is used for ghost element id etc */
}qid_ext_t;


// /**
//  * @brief defines the interaction of two sub set two field, i.e. a place in memory
//  * 
//  */
// typedef struct dom2dom_t {
//     // lid_t shift[3];
//     lid_t stride; /**< target memory stride */
//     lid_t gs;     /**< target ghost size */
//     lid_t start[3];
//     lid_t range[3];
//     // lid_t src_stride; /**< src memory stride */
//     // lid_t src_gs;     /**< src ghost size */
//     // lid_t src_start[3];
//     // lid_t src_range[3];
// } dom_t;

#endif  // SRC_MURPHY_HPP_
