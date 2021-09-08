#include "gridcallback.hpp"

#include <p8est_bits.h>

#include <list>

#include "grid/grid.hpp"
#include "grid/gridblock.hpp"
#include "tools/patch.hpp"
#include "operator/setvalues.hpp"

using std::string;
using std::list;

/**
 * @brief initiate a new block and store its address in the p4est quad
 * 
 * @warning no memory allocation is done at this point!
 */
void cback_CreateBlock(p8est_iter_volume_info_t* info, void* user_data) {
    m_begin;
    //-------------------------------------------------------------------------
    p8est_t*              forest     = info->p4est;
    p8est_quadrant_t*     quad       = info->quad;
    p4est_topidx_t        which_tree = info->treeid;
    p8est_connectivity_t* connect    = forest->connectivity;

    // sanity checks
    m_assert(sizeof(GridBlock*) == forest->data_size, "cannot cast the pointer, this is baaaad");

    // create the new block and store it's address
    real_t xyz[3];
    p8est_qcoord_to_vertex(connect, which_tree, quad->x, quad->y, quad->z, xyz);
    m_assert(quad->level >= 0, "the level=%d must be >=0", quad->level);
    real_t len = p4est_QuadLen(quad->level);
    p4est_SetGridBlock(quad, new GridBlock(len, xyz, quad->level));
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief delete a create block associated to a p4est quad
 */
void cback_DestroyBlock(p8est_iter_volume_info_t* info, void* user_data) {
    m_begin;
    //-------------------------------------------------------------------------
    p8est_quadrant_t* quad  = info->quad;
    GridBlock*        block = p4est_GetGridBlock(quad);
    m_assert(block != nullptr, "the block you are trying to free has already been free'ed");
    delete (block);
    p4est_SetGridBlock(quad, nullptr);
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief if within the level limits, always reply yes if p4est ask if we should refine the block
 */
int cback_Yes(p8est_t* forest, p4est_topidx_t which_tree, qdrt_t* quadrant) {
    m_begin;
    //-------------------------------------------------------------------------
    Grid* grid = reinterpret_cast<Grid*>(forest->user_pointer);

    // check the level limits
    bool is_refinable = grid->level_limit_min() <= quadrant->level && quadrant->level < grid->level_limit_max();

    if (is_refinable) {
        // add one block to the count, this drives the recursive adaptation
        grid->AddOneQuadToAdapt();
        return (true);
    } else {
        return (false);
    }
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief if within the level bounds, always reply yes if p4est ask if we should coarsen the block,
 */
int cback_Yes(p8est_t* forest, p4est_topidx_t which_tree, qdrt_t* quadrant[]) {
    m_begin;
    //-------------------------------------------------------------------------
    Grid* grid = reinterpret_cast<Grid*>(forest->user_pointer);

    
    bool is_coarsenable = true;
    for (iblock_t ib=0; ib< P8EST_CHILDREN; ++ib){
        is_coarsenable = is_coarsenable && (grid->level_limit_min() < quadrant[ib]->level) && (quadrant[ib]->level <= grid->level_limit_max());
    }

    if (is_coarsenable) {
        // add one block to the count, this drives the recursive adaptation
        grid->AddOneQuadToAdapt();
        return (true);
    } else {
        return (false);
    }
    //-------------------------------------------------------------------------
    m_end;
}

// /**
//  * @brief replies yes if we do not match one of the patch criterion
//  */
// int cback_Patch(p8est_t* forest, p4est_topidx_t which_tree, qdrt_t* quadrant) {
//     m_assert(false,"should go to trash");
//     // //-------------------------------------------------------------------------
//     // // retreive the patch list and the current block
//     // Grid*        grid    = reinterpret_cast<Grid*>(forest->user_pointer);
//     // list<Patch>* patches = reinterpret_cast<list<Patch>*>(grid->cback_criterion_ptr());
//     // GridBlock*   block   = p4est_GetGridBlock(quadrant);
//     // m_assert(block->level() == quadrant->level, "the two levels must match");
//     // m_assert(block->level() >= 0, "the level=%d must be >=0", block->level());

//     // // get the origin, the length and check if we are inside the patch
//     // const real_t* xyz = block->xyz();
//     // real_t        len = p4est_QuadLen(block->level());

//     // for (auto iter = patches->begin(); iter != patches->end(); iter++) {
//     //     Patch* patch = &(*iter);

//     //     // check that the patch's level is within the possible bounds
//     //     m_assert(grid->level_limit_min() <= patch->level(), "The patch level = %d must be >= min level = %d", patch->level(), grid->level_limit_min());
//     //     m_assert(patch->level() <= grid->level_limit_max(), "The patch level = %d must be <= max level = %d", patch->level(), grid->level_limit_max());

//     //     // if we already have the correct level or a higher one, we skip the patch
//     //     if (block->level() >= patch->level()) {
//     //         continue;
//     //     }
//     //     // if not, we have a coarser block and we might want to refine if the location matches
//     //     bool refine = true;
//     //     for (lda_t id = 0; id < 3; id++) {
//     //         // we have to satisfy both the our max > min and the min < our max
//     //         refine = refine &&
//     //                  (block->xyz(id) < (patch->origin(id) + patch->length(id))) &&
//     //                  (patch->origin(id) < (block->xyz(id) + len));
//     //     }
//     //     if (refine) {
//     //         grid->AddOneQuadToAdapt();
//     //         //m_log("adapt please!");
//     //         block->status_level(M_ADAPT_FINER);
//     //         return true;
//     //     }
//     // }

//     // // m_log("never mind");
//     // return false;
//     // //-------------------------------------------------------------------------
// }

// /**
//  * @brief replies yes if we do not match one of the patch criterion
//  */
// int cback_Patch(p8est_t* forest, p4est_topidx_t which_tree, qdrt_t* quadrant[]) {
//     m_assert(false,"should go to trash");
//     // //-------------------------------------------------------------------------
//     // // retreive the patch list and the current block
//     // Grid*        grid    = reinterpret_cast<Grid*>(forest->user_pointer);
//     // list<Patch>* patches = reinterpret_cast<list<Patch>*>(grid->cback_criterion_ptr());

//     // // check every block, if one child needs to be coarsen, we return true for everybody
//     // bool coarsen_all = false;
//     // for (sid_t ib = 0; ib < P8EST_CHILDREN; ib++) {
//     //     GridBlock* block = p4est_GetGridBlock(quadrant[ib]);
//     //     m_assert(block->level() == quadrant[ib]->level, "the two levels must match");
//     //     m_assert(block->level() >= 0, "the level=%d must be >=0", block->level());

//     //     // get the origin, the length and check if we are inside the patch
//     //     const real_t* xyz = block->xyz();
//     //     real_t        len = p4est_QuadLen(block->level());

//     //     for (auto iter = patches->begin(); iter != patches->end(); iter++) {
//     //         Patch* patch = &(*iter);

//     //         // check that the patch's level is within the possible bounds
//     //         m_assert(grid->level_limit_min() <= patch->level(), "The patch level = %d must be >= min level = %d", patch->level(), grid->level_limit_min());
//     //         m_assert(patch->level() <= grid->level_limit_max(), "The patch level = %d must be <= max level = %d", patch->level(), grid->level_limit_max());

//     //         m_log("comparing level: %d vs %d", block->level(), patch->level());
//     //         // if we already have the correct level or a lower one, we skip the patch
//     //         if (block->level() <= patch->level()) {
//     //             continue;
//     //         }
//     //         // if not, we have a finer block and we might want to coarsen if the location matches
//     //         bool coarsen = true;
//     //         for (lda_t id = 0; id < 3; id++) {
//     //             // we have to satisfy both the our max > min and the min < our max
//     //             coarsen = coarsen &&
//     //                       (block->xyz(id) < (patch->origin(id) + patch->length(id))) &&
//     //                       (patch->origin(id) < (block->xyz(id) + len));
//     //         }
//     //         // if only one needs to be coarsen, we coarsen all of them
//     //         if (coarsen) {
//     //             grid->AddOneQuadToAdapt();
//     //             coarsen_all = true;
//     //             break;
//     //         }
//     //     }
//     // }
//     // if (coarsen_all) {
//     //     for (sid_t ib = 0; ib < P8EST_CHILDREN; ib++) {
//     //         GridBlock* block = p4est_GetGridBlock(quadrant[ib]);
//     //         block->status_level(M_ADAPT_COARSER);
//     //     }
//     // }
//     // return coarsen_all;
//     // //-------------------------------------------------------------------------
// }

/**
 * @brief refine if the GridBlock::status_lvl_ of the associated block is +1
 */
int cback_StatusCheck(p8est_t* forest, p4est_topidx_t which_tree, qdrt_t* quadrant) {
    m_begin;
    //-------------------------------------------------------------------------
    Grid*      grid  = reinterpret_cast<Grid*>(forest->user_pointer);
    GridBlock* block = p4est_GetGridBlock(quadrant);

    const StatusAdapt current_status = block->status_level();

    // I allow refinement if we must and that the level is okay to do so (i.e. in the grid bounds)
    bool refine = current_status == M_ADAPT_FINER &&
                  grid->level_limit_min() <= block->level() &&
                  block->level() < grid->level_limit_max();

    if (refine) {
        grid->AddOneQuadToAdapt();
    } else {
        // update the status if we will not refine
        StatusAdapt new_status = (current_status == M_ADAPT_FINER) ? (M_ADAPT_SAME) : (current_status);
        block->status_level(new_status);
    }
    return refine;
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief coarsen if the GridBlock::status_lvl_ of every block is -1 and the level of each block is > level_limit_min
 * 
 * If one of the block needs to be refined or does not need to be coarsened, we cannot coarsen (safety first)
 */
int cback_StatusCheck(p8est_t* forest, p4est_topidx_t which_tree, qdrt_t* quadrant[]) {
    //-------------------------------------------------------------------------
    Grid* grid = reinterpret_cast<Grid*>(forest->user_pointer);

    // I allow coarsening if ALL the j8 children blocks are OK with it
    bool coarsen = true;
    // for each of the children, check if any of them prevent the coarsening
    for (short_t id = 0; id < P8EST_CHILDREN; id++) {
        GridBlock*        block          = p4est_GetGridBlock(quadrant[id]);
        const StatusAdapt current_status = block->status_level();
        // if we cannot coarsen for one of the block, return false
        coarsen = coarsen &&
                  (current_status == M_ADAPT_COARSER) &&
                  (block->level() > grid->level_limit_min());
    }
    // if I can coarsen the whole group, register them
    if (coarsen) {
        grid->AddQuadToAdapt(P8EST_CHILDREN);
    } else {
        // if not, make sure that the block that wanted to get coarsened have their tag changed
        for (short_t id = 0; id < P8EST_CHILDREN; id++) {
            GridBlock*        block          = p4est_GetGridBlock(quadrant[id]);
            const StatusAdapt current_status = block->status_level();
            const StatusAdapt new_status     = (current_status == M_ADAPT_COARSER) ? (M_ADAPT_SAME) : (current_status);
            block->status_level(new_status);
        }
    }
    return coarsen;
    //-------------------------------------------------------------------------
}

/**
 * @brief Update the depedency list of the incoming block(s) by storing the ID of the outgoing one(s).
 * 
 * The depedency list will be interpolated once the grid structure is fixed, here we just register the link
 * 
 * @param forest 
 * @param which_tree 
 * @param num_outgoing 
 * @param outgoing 
 * @param num_incoming 
 * @param incoming 
 */
void cback_UpdateDependency(p8est_t* forest, p4est_topidx_t which_tree, int num_outgoing, qdrt_t* outgoing[], int num_incoming, qdrt_t* incoming[]) {
    m_begin;
    m_assert(num_incoming == 1 || num_outgoing == 1, "we have either to compress or to refine");
    m_assert(num_incoming == P8EST_CHILDREN || num_outgoing == P8EST_CHILDREN, "the number of replacing blocks has to be the number of children");
    m_assert(forest->user_pointer != nullptr, "we need the grid in this function");
    //-------------------------------------------------------------------------
    // retrieve the grid from the forest user-data pointer
    Grid* grid = reinterpret_cast<Grid*>(forest->user_pointer);
    m_assert(grid->interp() != nullptr, "a Grid Wavelet is needed");
    // m_assert(!grid->recursive_adapt(), "the dependency update does not support recursive adaptation as we eventually need the GP values");

    // get needed grid info
    p8est_connectivity_t* connect = forest->connectivity;

    // count the number of blocks that have already been registered as dependencies
    iblock_t n_active_total = 0;
    for (iblock_t iout = 0; iout < num_outgoing; iout++) {
        GridBlock* block_out = p4est_GetGridBlock(outgoing[iout]);
        n_active_total += block_out->n_dependency_active();
    }
    m_assert(n_active_total == 0 || n_active_total == P8EST_CHILDREN, "the number of existing dependences (=%d) must be 0 or P8EST_CHILDREN", n_active_total);

    // create the incoming block if we have no existing blocks
    // if we already have registered depts, the GridBlocks already exist, so no need to do it
    for (sid_t iin = 0; iin < (num_incoming * (n_active_total == 0)); ++iin) {
        qdrt_t* quad = incoming[iin];

        // get block informations and create it
        real_t xyz[3];
        p8est_qcoord_to_vertex(connect, which_tree, quad->x, quad->y, quad->z, xyz);
        m_assert(quad->level >= 0, "the level=%d must be >=0", quad->level);
        real_t     len      = p4est_QuadLen(quad->level);
        GridBlock* block_in = new GridBlock(len, xyz, quad->level);
 
        // store the block in the quad
        p4est_SetGridBlock(quad, block_in);
    }

    // for each outgoing block
    for (sid_t iout = 0; iout < num_outgoing; iout++) {
        GridBlock* block_out = p4est_GetGridBlock(outgoing[iout]);
        // // check the status of the outgoing block, if from the balacing the status can be "do nothing"!
        // m_assert(!(num_incoming == 1 && (block_out->status_level() != M_ADAPT_COARSER || block_out->status_level() != M_ADAPT_SAME)), "the tag must be coarser instead of %d", block_out->status_level());
        // m_assert(!(num_incoming == P8EST_CHILDREN && (block_out->status_level() != M_ADAPT_FINER || block_out->status_level() != M_ADAPT_SAME)), "(num in = %d) an dthe tag must be finer instead of %d", num_incoming, block_out->status_level());
        // if called via the balancing it can be that the status doesn't match (of the outgoing block)

        // count the number of active dependency blocks that are already there
        sid_t n_active = block_out->n_dependency_active();
        m_assert(n_active == 0 || n_active == 1 || n_active == P8EST_CHILDREN, "the number of active dependency should always be 0 or %d, now: %d", P8EST_CHILDREN, n_active);

        if (n_active == 0) { /* most common: if we had no active dependencies, allocate new blocks, lock them and add them */
                             // we cannot have an old status or a none status here
                             // the status might be SAME because of the 2:1 imposition done by p4est
#ifndef NDEBUG
            // I can only refine an already refined  coarsen an already refined
            bool is_status_ok = (block_out->status_level() != M_ADAPT_NONE) &&
                                ((block_out->status_level() >= M_ADAPT_SAME) || !(num_incoming = 1 && block_out->status_level() == M_ADAPT_NEW_FINE));
            m_assert(is_status_ok, "the leaving block must have valid status: %d with incomming = %d", block_out->status_level(), incoming);
#endif

            // globaly we created the block ourselves
            m_assert(n_active_total == 0, "we must have created the associated block");
            for (sid_t iin = 0; iin < num_incoming; ++iin) {
                GridBlock* block_in = p4est_GetGridBlock(incoming[iin]);
                m_assert(block_in != nullptr, "block is null, ohoh");

                // assign the status
                StatusAdapt status_in = (num_incoming == 1) ? M_ADAPT_NEW_COARSE : M_ADAPT_NEW_FINE;
                block_in->status_level(status_in);

                // register the leaving block in the new one, using it's child id
                int childid_out = (num_outgoing == 1) ? 0 : p4est_GetChildID(block_out->xyz(), block_out->level());
                block_in->PushDependency(childid_out, block_out);

                // register the new block in the leaving one, using it's child id
                int childid_in = (num_incoming == 1) ? 0 : p4est_GetChildID(block_in->xyz(), block_in->level());
                block_out->PushDependency(childid_in, block_in);
                // assign the status
                StatusAdapt status_out = (num_incoming == 1) ? M_ADAPT_COARSER : M_ADAPT_FINER;
                block_out->status_level(status_out);
            }
        } else if (n_active == 1) { /* in case of 2:1 balacing: we want to coarsen blocks that have been asked for refinement */
            m_assert(n_active_total == P8EST_CHILDREN, "we must found a total of P8EST_CHILDREN blocks");
            m_assert(num_incoming == 1 && num_outgoing == P8EST_CHILDREN, "we cannot refine a block that has already been targeted for refinement!");

            // recover the old coarse block and remove the dependency from the block
            // this might be done several times, no problem
            GridBlock* parent = block_out->PopDependency(0);
            p4est_SetGridBlock(incoming[0], parent);

            // clear my name from my parent's list
            m_assert(p4est_GetChildID(block_out->xyz(), block_out->level()) == p8est_quadrant_child_id(outgoing[iout]), "ouuups");
            int childid = p4est_GetChildID(block_out->xyz(), block_out->level());
            parent->PopDependency(childid);

            // check the status, the parent must have been asked to refine and reset it
            m_assert(parent->status_level() == M_ADAPT_FINER, "the parent must not change its status (now %d)", parent->status_level());
            parent->status_level(M_ADAPT_SAME);

            // delete the created block, no need to update the status
            delete block_out;

        } else if (n_active == P8EST_CHILDREN) { /* we want to refine blocks that have been asked for coarsening */
            m_assert(n_active_total == P8EST_CHILDREN, "we must found a total of P8EST_CHILDREN blocks");
            m_assert(num_incoming == P8EST_CHILDREN && num_outgoing == 1, "we cannot refine a block that has already been targeted for refinement!");

            for (sid_t iin = 0; iin < num_incoming; iin++) {
                // get the correct child id and pop the existing block out
                int        childid = p8est_quadrant_child_id(incoming[iin]);
                GridBlock* child   = block_out->PopDependency(childid);

                // check the status (the child must have been asked to coarsen) and reset it
                m_assert(child->status_level() == M_ADAPT_COARSER, "the parent must not change its status (now %d)", child->status_level());
                child->status_level(M_ADAPT_SAME);

                // store the child adress as the new block
                p4est_SetGridBlock(incoming[iin], child);
                // clear y name from the child list
                child->PopDependency(0);
            }
            delete block_out;
        }
    }
    //-------------------------------------------------------------------------
    m_end;
};

/**
 * @brief allocate the entering children or allocate the entering parent, no interpolation is performed
 * 
 * @warning the existing fields are reset to 0.0
 * 
 * @param forest 
 * @param which_tree 
 * @param num_outgoing 
 * @param outgoing 
 * @param num_incoming 
 * @param incoming 
 */
void cback_AllocateOnly(p8est_t* forest, p4est_topidx_t which_tree, int num_outgoing, qdrt_t* outgoing[], int num_incoming, qdrt_t* incoming[]) {
    m_begin;
    m_assert(num_incoming == 1 || num_outgoing == 1, "we have either to compress or to refine");
    m_assert(num_incoming == P8EST_CHILDREN || num_outgoing == P8EST_CHILDREN, "the number of replacing blocks has to be the number of children");
    m_assert(forest->user_pointer != nullptr, "we need the grid in this function");
    //-------------------------------------------------------------------------
    // retrieve the grid from the forest user-data pointer
    Grid* grid = reinterpret_cast<Grid*>(forest->user_pointer);
    m_assert(grid->interp() != nullptr, "a Grid Wavelet is needed");

    // get needed grid info
    auto                  f_start = grid->FieldBegin();
    auto                  f_end   = grid->FieldEnd();
    p8est_connectivity_t* connect = forest->connectivity;

    // we need to know if we are re-assigning block that have already been flaged for refinement/coarsening
    iblock_t n_reassign = 0;
    for (iblock_t iout = 0; iout < num_outgoing; iout++) {
        GridBlock* block_out = p4est_GetGridBlock(outgoing[iout]);
        m_assert(block_out->status_level() != M_ADAPT_NONE, "the block cannot have a none status");
        n_reassign += (block_out->status_level() == M_ADAPT_NEW_COARSE) || (block_out->status_level() == M_ADAPT_NEW_FINE);
    }
    m_assert(n_reassign == 0 || n_reassign == num_outgoing, "the number of reassigning (=%d) must be 0 or %d", n_reassign, num_outgoing);

    // allocate the incomming blocks
    for (int id = 0; id < num_incoming; id++) {
        qdrt_t* quad = incoming[id];
        // get block informations and create it
        real_t xyz[3];
        p8est_qcoord_to_vertex(connect, which_tree, quad->x, quad->y, quad->z, xyz);
        m_assert(quad->level >= 0, "the level=%d must be >=0", quad->level);
        real_t     len   = p4est_QuadLen(quad->level);
        GridBlock* block = new GridBlock(len, xyz, quad->level);
        // store the block
        p4est_SetGridBlock(quad, block);
        // for every field, we allocate the memory
        for (auto fid = f_start; fid != f_end; fid++) {
            // allocate the new field
            block->AddField(fid->second);
        }

        // assign the status
        StatusAdapt status = (n_reassign == num_outgoing) ? (M_ADAPT_SAME) : ((num_incoming == 1) ? M_ADAPT_NEW_COARSE : M_ADAPT_NEW_FINE);
        block->status_level(status);
    }
    // deallocate the leaving blocks
    for (int id = 0; id < num_outgoing; id++) {
        qdrt_t* quad = outgoing[id];
        // GridBlock* block = *(reinterpret_cast<GridBlock**>(quad->p.user_data));
        GridBlock* block = p4est_GetGridBlock(quad);
        // delete the block, the fields are destroyed in the destructor
        delete (block);
    }
    //-------------------------------------------------------------------------
    m_end;
}

void cback_ValueFill(p8est_t* forest, p4est_topidx_t which_tree, int num_outgoing, qdrt_t* outgoing[], int num_incoming, qdrt_t* incoming[]) {
    m_begin;
    m_assert(num_incoming == 1 || num_outgoing == 1, "we have either to compress or to refine");
    m_assert(num_incoming == P8EST_CHILDREN || num_outgoing == P8EST_CHILDREN, "the number of replacing blocks has to be the number of children");
    m_assert(forest->user_pointer != nullptr, "we need the grid in this function");
    //-------------------------------------------------------------------------
    // retrieve the grid from the forest user-data pointer
    Grid*           grid  = static_cast<Grid*>(forest->user_pointer);
    Field*          field = static_cast<Field*>(grid->cback_criterion_ptr());
    SetValue* expr  = static_cast<SetValue*>(grid->cback_interpolate_ptr());

    m_assert(grid->interp() != nullptr, "a Grid Wavelet is needed");
    // m_assert(expr->DoGhost(), "the SetValue object must set the ghost values");
    m_assert(field != nullptr, "the field to fill shouldn't be nullptr");
    m_assert(grid->IsAField(field), "the field to fill must exist");

    // get needed grid info
    auto                  f_start = grid->FieldBegin();
    auto                  f_end   = grid->FieldEnd();
    p8est_connectivity_t* connect = forest->connectivity;

    // we need to know if we are re-assigning block that have already been flaged for refinement/coarsening
    iblock_t n_reassign = 0;
    for (iblock_t iout = 0; iout < num_outgoing; iout++) {
        GridBlock* block_out = p4est_GetGridBlock(outgoing[iout]);
        m_assert(block_out->status_level() != M_ADAPT_NONE, "the block cannot have a none status");
        n_reassign += (block_out->status_level() == M_ADAPT_NEW_COARSE) || (block_out->status_level() == M_ADAPT_NEW_FINE);
    }
    m_assert(n_reassign == 0 || n_reassign == num_outgoing, "the number of reassigning (=%d) must be 0 or %d", n_reassign, num_outgoing);

    // allocate the incomming blocks
    for (int id = 0; id < num_incoming; id++) {
        qdrt_t* quad = incoming[id];
        // get block informations and create it
        real_t xyz[3];
        p8est_qcoord_to_vertex(connect, which_tree, quad->x, quad->y, quad->z, xyz);
        m_assert(quad->level >= 0, "the level=%d must be >=0", quad->level);
        real_t     len   = p4est_QuadLen(quad->level);
        GridBlock* block = new GridBlock(len, xyz, quad->level);
        // store the block
        p4est_SetGridBlock(quad, block);
        // for every field, we allocate the memory
        for (auto fid = f_start; fid != f_end; fid++) {
            // allocate the new field
            block->AddField(fid->second);
        }
        // fill the new block with the analytical value
        expr->FillGridBlock(nullptr, block, field);

        // assign the status
        StatusAdapt status = (n_reassign == num_outgoing) ? (M_ADAPT_SAME) : ((num_incoming == 1) ? M_ADAPT_NEW_COARSE : M_ADAPT_NEW_FINE);
        block->status_level(status);
    }

    // deallocate the leaving blocks
    for (int id = 0; id < num_outgoing; id++) {
        qdrt_t*    quad  = outgoing[id];
        GridBlock* block = p4est_GetGridBlock(quad);
        // delete the block, the fields are destroyed in the destructor
        delete (block);
    }
    //-------------------------------------------------------------------------
    m_end;
}

// void cback_MGCreateFamilly(p8est_t* forest, p4est_topidx_t which_tree, int num_outgoing, qdrt_t* outgoing[], int num_incoming, qdrt_t* incoming[]) {
//     m_begin;
//     m_assert(num_outgoing == P8EST_CHILDREN, "this function is only called when doing the refinement");
//     m_assert(num_incoming == 1, "this function is only called when doing the refinement");
//     //-------------------------------------------------------------------------
//     p8est_connectivity_t* connect = forest->connectivity;
//     // retrieve the multigrid
//     Multigrid* grid = reinterpret_cast<Multigrid*>(forest->user_pointer);

//     // get the new block informations and create it
//     qdrt_t* quad = incoming[0];
//     real_t xyz[3];
//     p8est_qcoord_to_vertex(connect, which_tree, quad->x, quad->y, quad->z, xyz);
//     real_t     len      = p4est_QuadLen(quad->level);
//     GridBlock* parent = new GridBlock(len, xyz, quad->level);
//     // allocate the correct fields
//     parent->AddFields(grid->map_fields());
//     // store the new block
//     *(reinterpret_cast<GridBlock**>(quad->p.user_data)) = parent;

//     // get the leaving blocks
//     GridBlock* children[P8EST_CHILDREN];
//     for(sid_t ic=0; ic<P8EST_CHILDREN; ic++){
//         children[ic] = *(reinterpret_cast<GridBlock**>(outgoing[ic]->p.user_data));
//     }

//     // bind the family together
//     MGFamily* family = grid->curr_family();
//     family->AddMembers(parent,children);
//     //-------------------------------------------------------------------------
//     m_end;
// }

// /**
//  * @brief always reply yes if p4est ask if we should coarsen the block
//  */
// int cback_Level(p8est_t* forest, p4est_topidx_t which_tree, qdrt_t* quadrant[]) {
//     m_begin;
//     //-------------------------------------------------------------------------
//     Multigrid* grid = reinterpret_cast<Multigrid*>(forest->user_pointer);
//     sid_t target_level = grid->curr_level();
//     return (quadrant[0]->level > target_level);
//     //-------------------------------------------------------------------------
//     m_end;
// }