#ifndef SRC_MURPHY_HPP_
#define SRC_MURPHY_HPP_

#include <p8est.h>

#include <map>
#include <string>

#include "defs.hpp"

typedef int32_t lid_t;
typedef int8_t  sid_t;

/**
 * @name typename for indexes
 * @{
 */
typedef int8_t         iface_t;  //!< face index type (0->26)
typedef int8_t         lda_t;    //!< leading dimension array type
typedef int8_t         level_t;  //!< data type for a level
typedef int            rank_t;   //!< rank data type
typedef p4est_locidx_t iblock_t;  //!< datatype to index blocks, ghosts and quadrants
/**@}*/

/**
 * @name typename for data pointers
 * @{
 */
typedef double real_t;              //!< data subtype
typedef real_t* __restrict real_p;  //!< pointer type = root of the memory allocation
typedef real_t* __restrict data_p;  //!< data pointer type = root of the data, i.e. point (0,0,0)
/**@}*/

typedef std::map<std::string, real_p> datamap_t;

typedef p8est_quadrant_t qdrt_t;

/**
 * @brief id of a quadrant
 * 
 */
typedef struct qid_t
{
    p4est_locidx_t cid; //!< cummulative_id, does not depends on the tree id
    p4est_locidx_t qid; //!< quadrant id = depends on the tree
    p4est_topidx_t tid; //!< tree id
} qid_t;

/**
 * @brief defines the different supported boundary conditions
 * 
 * The values matches the FLUPS values
 * 
 */
typedef enum bctype_t
{
    M_BC_NONE,     //!< no boundary condition is given
    M_BC_EVEN,     //!< EVEN condition: an EVEN bondary condition with respect to a given flux at the interface
    M_BC_ODD,      //!< ODD condition: imposes a zero value wrt the interface
    M_BC_ZERO,     //!< set 0 outside the computational domain
    M_BC_EXTRAP_3, //!< extrapolate outside the computational domain with polynomial x^2
    M_BC_EXTRAP_4, //!< extrapolate outside the computational domain with polynomial x^3
    M_BC_EXTRAP_5  //!< extrapolate outside the computational domain with polynomial x^4
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
