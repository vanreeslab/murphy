#include "gridcallback.hpp"

#include <p8est_bits.h>
#include <list>

#include "grid/cartblock.hpp"
#include "grid/grid.hpp"
#include "grid/gridblock.hpp"
#include "tools/toolsp4est.hpp"
#include "tools/toolsblocktypes.hpp"

using std::list;
using std::string;

/**
 * @brief Returns the associated cback_CreateBlock function for p4est
 */
p8est_iter_volume_t get_cback_CreateBlock(BlockDataType block_type) {
    m_assert(block_type != M_NULLTYPE, "cannot provide size of unspecified block type.");
    if (block_type == M_CARTBLOCK) {
        return cback_CreateBlock<CartBlock>;
    }
    if (block_type == M_GRIDBLOCK) {
        return cback_CreateBlock<GridBlock>;
    }
    // if you get here, there is an error
    m_assert(false, "BlockDataType value not recognized.");
    return (p8est_iter_volume_t)(nullptr);
}

/**
 * @brief Returns the associated cback_DestroyBlock function for p4est
 */
p8est_iter_volume_t get_cback_DestroyBlock(BlockDataType block_type) {
    m_assert(block_type != M_NULLTYPE, "cannot provide size of unspecified block type.");
    if (block_type == M_CARTBLOCK) {
        return cback_DestroyBlock<CartBlock>;
    }
    if (block_type == M_GRIDBLOCK) {
        return cback_DestroyBlock<GridBlock>;
    }
    // if you get here, there is an error
    m_assert(false, "BlockDataType value not recognized.");
    return (p8est_iter_volume_t)(nullptr);
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
        grid->AddOneQuadToRefine();
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
    for (iblock_t ib = 0; ib < P8EST_CHILDREN; ++ib) {
        is_coarsenable = is_coarsenable && (grid->level_limit_min() < quadrant[ib]->level) && (quadrant[ib]->level <= grid->level_limit_max());
    }

    if (is_coarsenable) {
        // add one block to the count, this drives the recursive adaptation
        grid->AddOneQuadToCoarsen();
        return (true);
    } else {
        return (false);
    }
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief refine given the GridBlock::status_lvl()
 */
int cback_StatusCheck(p8est_t* forest, p4est_topidx_t which_tree, qdrt_t* quadrant) {
    m_begin;
    //-------------------------------------------------------------------------
    Grid* grid  = reinterpret_cast<Grid*>(forest->user_pointer);
    auto  block = p4est_GetBlock<GridBlock>(quadrant);

    const StatusAdapt current_status = block->status_level();

    // I allow refinement if we must and that the level is okay to do so (i.e. in the grid bounds)
    bool refine = current_status == M_ADAPT_FINER &&
                  grid->level_limit_min() <= block->level() &&
                  block->level() < grid->level_limit_max();

    if (refine) {
        grid->AddOneQuadToRefine();
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
    m_assert(IsCompatibleBlockType(grid->block_type(), TypeToEnum<GridBlock>()), "the grid must support a GridBlock");

    // I allow coarsening if ALL the j8 children blocks are OK with it
    bool coarsen = true;
    // for each of the children, check if any of them prevent the coarsening
    for (short_t id = 0; id < P8EST_CHILDREN; id++) {
        auto        block          = p4est_GetBlock<GridBlock>(quadrant[id]);
        const StatusAdapt current_status = block->status_level();
        // if we cannot coarsen for one of the block, return false
        coarsen = coarsen &&
                  (current_status == M_ADAPT_COARSER) &&
                  (block->level() > grid->level_limit_min());
    }
    // if I can coarsen the whole group, register them
    if (coarsen) {
        grid->AddOneQuadToCoarsen();
    } else {
        // if not, make sure that the block that wanted to get coarsened have their tag changed
        for (short_t id = 0; id < P8EST_CHILDREN; id++) {
            auto        block          = p4est_GetBlock<GridBlock>(quadrant[id]);
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
    m_assert(IsCompatibleBlockType(grid->block_type(),TypeToEnum<GridBlock>()),"The grid must support at least GridBlock to go through this routine");

    // get needed grid info
    p8est_connectivity_t* connect = forest->connectivity;

    // count the number of blocks that have already been registered as dependencies
    iblock_t n_active_total = 0;
    for (iblock_t iout = 0; iout < num_outgoing; iout++) {
        auto block_out = p4est_GetBlock<GridBlock>(outgoing[iout]);
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
        AllocateBlockForP4est(grid->block_type(), xyz, quad);
    }

    // for each outgoing block
    for (sid_t iout = 0; iout < num_outgoing; iout++) {
        auto block_out = p4est_GetBlock<GridBlock>(outgoing[iout]);

        // count the number of active dependency blocks that are already there
        sid_t n_active = block_out->n_dependency_active();
        m_assert(n_active == 0 || n_active == 1 || n_active == P8EST_CHILDREN, "the number of active dependency should always be 0 or %d, now: %d", P8EST_CHILDREN, n_active);

        if (n_active == 0) { /* most common: if we had no active dependencies, allocate new blocks, lock them and add them */
            // globaly we created the block ourselves
            m_assert(n_active_total == 0, "we must have created the associated block");
            for (sid_t iin = 0; iin < num_incoming; ++iin) {
                auto block_in = p4est_GetBlock<GridBlock>(incoming[iin]);
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
            p4est_SetBlock(incoming[0], parent);

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
                p4est_SetBlock(incoming[iin], child);
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
    m_assert(IsCompatibleBlockType(grid->block_type(),TypeToEnum<GridBlock>()),"The grid must support at least GridBlock to go through this routine");

    // get needed grid info
    auto                  f_start = grid->FieldBegin();
    auto                  f_end   = grid->FieldEnd();
    p8est_connectivity_t* connect = forest->connectivity;

    // we need to know if we are re-assigning block that have already been flaged for refinement/coarsening
    iblock_t n_reassign = 0;
    for (iblock_t iout = 0; iout < num_outgoing; iout++) {
        auto block_out = p4est_GetBlock<GridBlock>(outgoing[iout]);
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
        AllocateBlockForP4est(grid->block_type(), xyz, quad);
        auto block = p4est_GetBlock<GridBlock>(quad);

        // for every field, we allocate the memory
        for (auto fid = f_start; fid != f_end; fid++) {
            // allocate the new field
            auto current_fid = fid->second;
            if (!current_fid->is_expr()) {
                block->AddField(current_fid);
            }
        }

        // assign the status
        StatusAdapt status = (n_reassign == num_outgoing) ? (M_ADAPT_SAME) : ((num_incoming == 1) ? M_ADAPT_NEW_COARSE : M_ADAPT_NEW_FINE);
        block->status_level(status);
    }
    // deallocate the leaving blocks
    for (int id = 0; id < num_outgoing; id++) {
        auto block = p4est_GetBlock<GridBlock>(outgoing[id]);
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
    Grid*     grid  = static_cast<Grid*>(forest->user_pointer);
    Field*    field = static_cast<Field*>(grid->cback_criterion_ptr());
    SetValue* expr  = static_cast<SetValue*>(grid->cback_interpolate_ptr());

    m_assert(grid->interp() != nullptr, "a Grid Wavelet is needed");
    m_assert(IsCompatibleBlockType(grid->block_type(),TypeToEnum<GridBlock>()),"The grid must support at least GridBlock to go through this routine");
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
        auto block_out = p4est_GetBlock<GridBlock>(outgoing[iout]);
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
        AllocateBlockForP4est(grid->block_type(), xyz, quad);
        auto block = p4est_GetBlock<GridBlock>(quad);
        // for every field, we allocate the memory
        for (auto fid = f_start; fid != f_end; fid++) {
            // allocate the new field
            auto current_fid = fid->second;
            if (!current_fid->is_expr()) {
                block->AddField(current_fid);
            }
        }
        // fill the new block with the analytical value
        expr->FillGridBlock(nullptr, block, field);

        // assign the status
        StatusAdapt status = (n_reassign == num_outgoing) ? (M_ADAPT_SAME) : ((num_incoming == 1) ? M_ADAPT_NEW_COARSE : M_ADAPT_NEW_FINE);
        block->status_level(status);
    }

    // deallocate the leaving blocks and set nullptr instead
    for (int id = 0; id < num_outgoing; id++) {
        qdrt_t*    quad  = outgoing[id];
        auto block = p4est_GetBlock<GridBlock>(quad);
        // delete the block, the fields are destroyed in the destructor
        delete (block);
        p4est_SetBlock<std::nullptr_t>(quad, nullptr);
    }
    //-------------------------------------------------------------------------
    m_end;
}