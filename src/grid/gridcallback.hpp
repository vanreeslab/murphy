#ifndef SRC_GRIDCALLBACK_HPP_
#define SRC_GRIDCALLBACK_HPP_

#include <p8est.h>
#include <p8est_extended.h>
#include "core/macros.hpp"
#include "core/types.hpp"
#include "tools/toolsp4est.hpp"

using cback_refine_criterion_t = int (*)(p8est_t *forest, p4est_topidx_t which_tree, qdrt_t *quadrant);
using cback_coarsen_citerion_t = int (*)(p8est_t *forest, p4est_topidx_t which_tree, qdrt_t *quadrant[]);
using cback_interpolate_t      = void (*)(p8est_t *forest, p4est_topidx_t which_tree, int num_outgoing, qdrt_t *outgoing[], int num_incoming, qdrt_t *incoming[]);


/**
 * @name create and destroy blocks
 * 
 * @{
 */

/**
 * @brief initiate a new block and store its address in the p4est quad
 * 
 * @warning no memory allocation is done at this point!
 */

template <typename BlockType>
void cback_CreateBlock(p8est_iter_volume_info_t *info, void *user_data) {
    m_begin;
    //-------------------------------------------------------------------------
    p8est_t *             forest     = info->p4est;
    p8est_quadrant_t *    quad       = info->quad;
    p4est_topidx_t        which_tree = info->treeid;
    p8est_connectivity_t *connect    = forest->connectivity;
    // sanity checks
    m_assert(sizeof(BlockType *) == forest->data_size, "cannot cast the pointer, this is baaaad");

    // create the new block and store it's address
    real_t xyz[3];
    p8est_qcoord_to_vertex(connect, which_tree, quad->x, quad->y, quad->z, xyz);
    m_assert(quad->level >= 0, "the level=%d must be >=0", quad->level);
    real_t     len   = p4est_QuadLen(quad->level);
    BlockType *block = new BlockType(len, xyz, quad->level);
    p4est_SetBlock(quad, block);
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief delete a create block associated to a p4est quad
 */
template <typename BlockType>
void cback_DestroyBlock(p8est_iter_volume_info_t *info, void *user_data) {
    m_begin;
    //-------------------------------------------------------------------------
    // get the quad
    p8est_quadrant_t *quad = info->quad;

    // retrieve the qad from the block and delete everything
    auto block = p4est_GetBlock<BlockType>(quad);
    m_assert(block != nullptr, "the block you are trying to free has already been free'ed");
    delete (block);
    p4est_SetBlock<std::nullptr_t>(quad, nullptr);
    //-------------------------------------------------------------------------
    m_end;
}

p8est_iter_volume_t get_cback_CreateBlock(BlockDataType block_type);
p8est_iter_volume_t get_cback_DestroyBlock(BlockDataType block_type);
/**@}*/

/**
 * @name refine and coarsen
 * 
 * @{
 */
// yes all the time
int cback_Yes(p8est_t *forest, p4est_topidx_t which_tree, qdrt_t *quadrant);
int cback_Yes(p8est_t *forest, p4est_topidx_t which_tree, qdrt_t *quadrant[]);
// check the status form the block
int cback_StatusCheck(p8est_t *forest, p4est_topidx_t which_tree, qdrt_t *quadrant);
int cback_StatusCheck(p8est_t *forest, p4est_topidx_t which_tree, qdrt_t *quadrant[]);

// interpolation
void cback_AllocateOnly(p8est_t *forest, p4est_topidx_t which_tree, int num_outgoing, qdrt_t *outgoing[], int num_incoming, qdrt_t *incoming[]);
void cback_ValueFill(p8est_t *forest, p4est_topidx_t which_tree, int num_outgoing, qdrt_t *outgoing[], int num_incoming, qdrt_t *incoming[]);
void cback_UpdateDependency(p8est_t *forest, p4est_topidx_t which_tree, int num_outgoing, qdrt_t *outgoing[], int num_incoming, qdrt_t *incoming[]);
/**@}*/

#endif  // SRC_GRIDCALLBACK_HPP_
