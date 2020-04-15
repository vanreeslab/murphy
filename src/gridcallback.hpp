#ifndef SRC_GRID_CALLBACK_HPP_
#define SRC_GRID_CALLBACK_HPP_

#include "murphy.hpp"
#include <p8est.h>
#include <p8est_extended.h>

/**
 * @name create and destroy blocks
 * 
 * @{
 */
void cback_CreateBlock(p8est_iter_volume_info_t* info, void* user_data);
void cback_DestroyBlock(p8est_iter_volume_info_t* info, void* user_data);
/**@}*/


/**
 * @name refine and coarsen
 * 
 * @{
 */
int cback_Yes(p8est_t *forest, p4est_topidx_t which_tree, qdrt_t *quadrant);
int cback_Yes(p8est_t *forest, p4est_topidx_t which_tree, qdrt_t *quadrant[]);
void cback_Interpolate(p8est_t *grid, p4est_topidx_t which_tree, int num_outgoing, qdrt_t *outgoing[], int num_incoming, qdrt_t *incoming[]);
/**@}*/



#endif