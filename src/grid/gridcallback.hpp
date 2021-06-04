#ifndef SRC_GRIDCALLBACK_HPP_
#define SRC_GRIDCALLBACK_HPP_

#include <p8est.h>
#include <p8est_extended.h>

#include "core/macros.hpp"
#include "core/types.hpp"

using cback_coarsen_citerion_t = int (*)(p8est_t *forest, p4est_topidx_t which_tree, qdrt_t *quadrant[]);
using cback_refine_criterion_t = int (*)(p8est_t *forest, p4est_topidx_t which_tree, qdrt_t *quadrant);
using cback_interpolate_t      = void (*)(p8est_t *forest, p4est_topidx_t which_tree, int num_outgoing, qdrt_t *outgoing[], int num_incoming, qdrt_t *incoming[]);

/**
 * @name create and destroy blocks
 * 
 * @{
 */
void cback_CreateBlock(p8est_iter_volume_info_t *info, void *user_data);
void cback_DestroyBlock(p8est_iter_volume_info_t *info, void *user_data);
/**@}*/

/**
 * @name refine and coarsen
 * 
 * @{
 */
// yes all the time
int cback_Yes(p8est_t *forest, p4est_topidx_t which_tree, qdrt_t *quadrant);
int cback_Yes(p8est_t *forest, p4est_topidx_t which_tree, qdrt_t *quadrant[]);
// yes if in the patch
// int cback_Patch(p8est_t *forest, p4est_topidx_t which_tree, qdrt_t *quadrant);
// int cback_Patch(p8est_t *forest, p4est_topidx_t which_tree, qdrt_t *quadrant[]);
// check the status form the block
int cback_StatusCheck(p8est_t *forest, p4est_topidx_t which_tree, qdrt_t *quadrant);
int cback_StatusCheck(p8est_t *forest, p4est_topidx_t which_tree, qdrt_t *quadrant[]);

// interpolation
// void cback_Interpolate(p8est_t *forest, p4est_topidx_t which_tree, int num_outgoing, qdrt_t *outgoing[], int num_incoming, qdrt_t *incoming[]);
void cback_AllocateOnly(p8est_t *forest, p4est_topidx_t which_tree, int num_outgoing, qdrt_t *outgoing[], int num_incoming, qdrt_t *incoming[]);
void cback_ValueFill(p8est_t *forest, p4est_topidx_t which_tree, int num_outgoing, qdrt_t *outgoing[], int num_incoming, qdrt_t *incoming[]);
void cback_UpdateDependency(p8est_t *forest, p4est_topidx_t which_tree, int num_outgoing, qdrt_t *outgoing[], int num_incoming, qdrt_t *incoming[]);
/**@}*/

/**
 * @name MultiGrid
 * 
 * @{
 */
// int  cback_Level(p8est_t *forest, p4est_topidx_t which_tree, qdrt_t *quadrant);
// int cback_Level(p8est_t *forest, p4est_topidx_t which_tree, qdrt_t *quadrant[]);
// void cback_ReplaceByMultiGridBlock(p8est_iter_volume_info_t *info, void *user_data);
// void cback_MGCreateFamilly(p8est_t *forest, p4est_topidx_t which_tree, int num_outgoing, qdrt_t *outgoing[], int num_incoming, qdrt_t *incoming[]);
/**@}*/

#endif  // SRC_GRIDCALLBACK_HPP_
