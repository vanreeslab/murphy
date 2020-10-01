#include "gridcallback.hpp"

#include <p8est_bits.h>

#include <list>

#include "grid.hpp"
#include "gridblock.hpp"
#include "patch.hpp"
#include "mgfamily.hpp"
#include "multigrid.hpp"
#include "setvalues.hpp"

using std::string;
using std::list;

/**
 * @brief initiate a new block and store its adress in the p4est quad
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
    // get the starting position
    real_t xyz[3];
    p8est_qcoord_to_vertex(connect, which_tree, quad->x, quad->y, quad->z, xyz);

    real_t len = m_quad_len(quad->level);
    // the user data points to the data defined by the grid = GridBlock*
    m_assert(sizeof(GridBlock*) == forest->data_size,"cannot cast the pointer, this is baaaad");
    *(reinterpret_cast<GridBlock**>(quad->p.user_data)) = new GridBlock(len, xyz, quad->level);
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
    GridBlock*        block = *(reinterpret_cast<GridBlock**>(quad->p.user_data));
    m_assert(block != nullptr,"the block you are trying to free has already been free'ed");
    delete (block);
    *(reinterpret_cast<GridBlock**>(quad->p.user_data)) = nullptr;
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief always reply yes if p4est ask if we should refine the block
 */
int cback_Yes(p8est_t* forest, p4est_topidx_t which_tree, qdrt_t* quadrant) {
    m_begin;
    //-------------------------------------------------------------------------
    return (true);
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief always reply yes if p4est ask if we should coarsen the block
 */
int cback_Yes(p8est_t* forest, p4est_topidx_t which_tree, qdrt_t* quadrant[]) {
    m_begin;
    //-------------------------------------------------------------------------
    return (true);
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief replies yes if we do not match one of the patch criterion
 */
int cback_Patch(p8est_t* forest, p4est_topidx_t which_tree, qdrt_t* quadrant) {
    //-------------------------------------------------------------------------
    // retreive the patch list and the current block
    Grid*        grid    = reinterpret_cast<Grid*>(forest->user_pointer);
    list<Patch>* patches = reinterpret_cast<list<Patch>*>(grid->cback_criterion_field());
    GridBlock*   block   = *(reinterpret_cast<GridBlock**>(quadrant->p.user_data));
    m_assert(block->level() == quadrant->level, "the two levels must match");

    // get the origin, the length and check if we are inside the patch
    const real_t* xyz = block->xyz();
    real_t len = m_quad_len(block->level());

    for(auto iter=patches->begin(); iter!= patches->end(); iter++){
        Patch* patch = &(*iter);
        
        // if we already have the correct level or a higher one, we skip the patch
        if(block->level() >= patch->level()){
            continue;
        }
        // if not, we have a coarser block and we might want to refine if the location matches
        bool refine = true;
        for (int id = 0; id < 3; id++) {
            // we have to satisfy both the our max > min and the min < our max
            refine = refine &&
                     (block->xyz(id) < (patch->origin(id) + patch->length(id))) &&
                     (patch->origin(id) < (block->xyz(id) + len));
        }
        // m_log("should be refined? %d: levels %d vs %d",refine,block->level(),patch->level());
        if(refine){
            return true;
        }
    }
    return false;
    //-------------------------------------------------------------------------
}

/**
 * @brief replies yes if we do not match one of the patch criterion
 */
int cback_Patch(p8est_t* forest, p4est_topidx_t which_tree, qdrt_t* quadrant[]) {
    //-------------------------------------------------------------------------
    // retreive the patch list and the current block
    Grid*        grid    = reinterpret_cast<Grid*>(forest->user_pointer);
    list<Patch>* patches = reinterpret_cast<list<Patch>*>(grid->cback_criterion_field());

    // check every block, if one child needs to be coarsen, we return true for everybody
    for (sid_t ib = 0; ib < P8EST_CHILDREN; ib++) {
        GridBlock* block = *(reinterpret_cast<GridBlock**>(quadrant[ib]->p.user_data));
        m_assert(block->level() == quadrant[ib]->level, "the two levels must match");

        // get the origin, the length and check if we are inside the patch
        const real_t* xyz = block->xyz();
        real_t        len = m_quad_len(block->level());

        for (auto iter = patches->begin(); iter != patches->end(); iter++) {
            Patch* patch = &(*iter);
            // if we already have the correct level or a lower one, we skip the patch
            if (block->level() <= patch->level()) {
                continue;
            }
            // if not, we have a finer block and we might want to coarsen if the location matches
            bool coarsen = false;
            for (sid_t id = 0; id < 3; id++) {
                // we have to satisfy both the our max > min and the min < our max
                coarsen = coarsen &&
                         (block->xyz(id) < (patch->origin(id) + patch->length(id))) &&
                         (patch->origin(id) < (block->xyz(id) + len));
            }
            if (coarsen) {
                return true;
            }
        }
    }
    return false;
    //-------------------------------------------------------------------------
}

/**
 * @brief refine a block if @ref Interpolator::Criterion() is bigger than @ref Grid::rtol()
 * 
 * @param forest 
 * @param which_tree 
 * @param quadrant 
 * @return int 
 */
int cback_Interpolator(p8est_t* forest, p4est_topidx_t which_tree, qdrt_t* quadrant) {
    m_begin;
    //-------------------------------------------------------------------------
    Grid*         grid   = reinterpret_cast<Grid*>(forest->user_pointer);
    Field*        fid    = reinterpret_cast<Field*>(grid->cback_criterion_field());
    GridBlock*    block  = *(reinterpret_cast<GridBlock**>(quadrant->p.user_data));
    Interpolator* interp = grid->interp();
    
    // if the block is locked, I cannot touch it anymore
    if (block->locked()) {
        m_log("block is locked, we do nothing!");
        return false;
    }
    m_log("let's check the refinement for the field %s",fid->name().c_str());

    bool refine = false;
    // if we are not locked, we are available for refinement
    // we refine the block if one of the field component needs it
    for (int ida = 0; ida < fid->lda(); ida++) {
        real_p data = block->data(fid, ida);
        real_t norm = interp->Criterion(block, data);
        // refine if the norm is bigger
        refine = (norm > grid->rtol());
        // refine the block if one dimension needs to be refined
        if (refine) {
            // lock the block indicate that it cannot be changed in res anymore
            block->lock();
            m_log("block %f %f %f: YES refine (%e > %e)", block->xyz(0), block->xyz(1), block->xyz(2), norm, grid->rtol());
            // return to refine
            return true;
        }
        else{
            m_log("block %f %f %f: NO refine (%e < %e)", block->xyz(0), block->xyz(1), block->xyz(2), norm, grid->rtol());
        }
    }
    // if we reached here, we do not need to refine
    return false;
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief reply that we should coarsen the group of 8 blocks if @ref Interpolator::Criterion() is lower than @ref Grid::ctol() for everyblock.
 * 
 * We do NOT coarsen if one of the block does not match the criterion
 * 
 * @param forest 
 * @param which_tree 
 * @param quadrant 
 * @return int 
 */
int cback_Interpolator(p8est_t* forest, p4est_topidx_t which_tree, qdrt_t* quadrant[]) {
    m_begin;
    //-------------------------------------------------------------------------
    Grid*         grid   = reinterpret_cast<Grid*>(forest->user_pointer);
    Field*        fid    = reinterpret_cast<Field*>(grid->cback_criterion_field());
    Interpolator* interp = grid->interp();

    // for each of the children
    bool coarsen = true;
    for (int id = 0; id < P8EST_CHILDREN; id++) {
        GridBlock* block = *(reinterpret_cast<GridBlock**>(quadrant[id]->p.user_data));

        // if one of the 8 block is locked, I cannot change it, neither the rest of the group
        if (block->locked()) {
            m_log("block is locked, we do nothing!");
            return false;
        }

        // check if I can coarsen
        for (int ida = 0; ida < fid->lda(); ida++) {
            real_p data = block->data(fid, ida);
            real_t norm = interp->Criterion(block, data);
            // coarsen if the norm is smaller than the tol
            coarsen = (norm < grid->ctol());
            m_log("coarsen? %e <? %e", norm, grid->ctol());

            // if I cannot coarsen, I can give up on the whole group, so return false
            if (!coarsen) {
                return false;
            }
        }
    }
    // if I arrived here, I can coarsen the whole group, so lock them and return true
    for (int id = 0; id < P8EST_CHILDREN; id++) {
        GridBlock* block = *(reinterpret_cast<GridBlock**>(quadrant[id]->p.user_data));
        block->lock();
    }
    return true;
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief Interpolate one block to his children or 8 children to their parent, using the Interpolator of the @ref Grid.
 * 
 * If one of the outgoing block(s) is locked, then I lock the incomming block(s)
 * 
 * @param forest 
 * @param which_tree 
 * @param num_outgoing 
 * @param outgoing 
 * @param num_incoming 
 * @param incoming 
 */
void cback_Interpolate(p8est_t* forest, p4est_topidx_t which_tree, int num_outgoing, qdrt_t* outgoing[], int num_incoming, qdrt_t* incoming[]) {
    m_begin;
    m_assert(num_incoming == 1 || num_outgoing == 1, "we have either to compress or to refine");
    m_assert(num_incoming == P8EST_CHILDREN || num_outgoing == P8EST_CHILDREN, "the number of replacing blocks has to be the number of children");
    m_assert(forest->user_pointer != nullptr, "we need the grid in this function");
    //-------------------------------------------------------------------------
    // retrieve the grid from the forest user-data pointer
    Grid* grid = reinterpret_cast<Grid*>(forest->user_pointer);
    m_assert(grid->interp() != nullptr, "a Grid interpolator is needed");
    m_assert(!grid->recursive_adapt(), "the wavelet refinement does not support recursive adaptation as no GP is filled");

    // get needed grid info
    auto                  f_start = grid->FieldBegin();
    auto                  f_end   = grid->FieldEnd();
    Interpolator*         interp  = grid->interp();
    p8est_connectivity_t* connect = forest->connectivity;

    // m_log("prof callback");
    m_profStart(grid->profiler(), "cback_interpolate");

    // m_log("interpolate callback");

    // allocate the incomming blocks
    for (int id = 0; id < num_incoming; id++) {
        qdrt_t* quad = incoming[id];
        // get block informations and create it
        real_t xyz[3];
        p8est_qcoord_to_vertex(connect, which_tree, quad->x, quad->y, quad->z, xyz);
        real_t     len   = m_quad_len(quad->level);
        GridBlock* block = new GridBlock(len, xyz, quad->level);
        // store the block
        *(reinterpret_cast<GridBlock**>(quad->p.user_data)) = block;
        // for every field, we allocate the memory
        for (auto fid = f_start; fid != f_end; fid++) {
            // allocate the new field
            block->AddField(fid->second);
        }
    }

    // if only one is leaving, it means we refine -> dlvl = -1
    // if only one is entering, it means we coarse -> dlvl = +1
    const sid_t dlvl = (num_outgoing == 1) ? -1 : 1;

    // do the required iterpolation
    for (sid_t iout = 0; iout < num_outgoing; iout++) {
        GridBlock* block_out = *(reinterpret_cast<GridBlock**>(outgoing[iout]->p.user_data));

        for (sid_t iin = 0; iin < num_incoming; iin++) {
            GridBlock* block_in = *(reinterpret_cast<GridBlock**>(incoming[iin]->p.user_data));

            // get the correct child id: if 1 is going out, the children are the incoming
            int childid = p8est_quadrant_child_id((num_outgoing == 1) ? incoming[iin] : outgoing[iout]);

            // if we refine
            if (num_outgoing == 1) {
                // check if the parent is locked and lock the child if this is the case
                // this may happen if the parent has been tagged for refinement
                if (block_out->locked()) {
                    block_in->lock();
                }

                // get the shift given the child id
                const lid_t shift[3] = {M_HN * ((childid % 2)), M_HN * ((childid % 4) / 2), M_HN * ((childid / 4))};

                // we create a subblock with the correct memory representation for the source = parent
                const lid_t src_start[3] = {shift[0] - M_GS, shift[1] - M_GS, shift[2] - M_GS};
                const lid_t src_end[3]   = {shift[0] + M_HN + M_GS, shift[1] + M_HN + M_GS, shift[2] + M_HN + M_GS};
                SubBlock    mem_src(M_GS, M_STRIDE, src_start, src_end);

                // for every field, we interpolate it
                for (auto fid = f_start; fid != f_end; fid++) {
                    Field* current_field = fid->second;
                    // interpolate for every dimension
                    for (sid_t ida = 0; ida < current_field->lda(); ida++) {
                        // get the pointers
                        interp->Interpolate(dlvl, shift, &mem_src, block_out->data(current_field, ida), block_in, block_in->data(current_field, ida));
                    }
                }
            } else if (num_incoming == 1) {
                // check if the child is locked and lock the parent if this is the case
                // this may happen if the parent has been tagged for refinement
                if (block_out->locked()) {
                    block_in->lock();
                }

                // get the shift
                const lid_t shift[3] = {-M_N * ((childid % 2)), -M_N * ((childid % 4) / 2), -M_N * ((childid / 4))};

                // we create a subblock with the correct memory representation for the target
                const lid_t trg_start[3] = {M_HN * ((childid % 2)), M_HN * ((childid % 4) / 2), M_HN * ((childid / 4))};
                const lid_t trg_end[3]   = {trg_start[0] + M_HN, trg_start[1] + M_HN, trg_start[2] + M_HN};
                SubBlock    mem_trg(M_GS, M_STRIDE, trg_start, trg_end);
                // and an extended source block
                const lid_t src_start[3] = {-M_GS, -M_GS, -M_GS};
                const lid_t src_end[3]   = {M_N + M_GS, M_N + M_GS, M_N + M_GS};
                SubBlock    mem_src(M_GS, M_STRIDE, src_start, src_end);

                // for every field, we interpolate it
                for (auto fid = f_start; fid != f_end; fid++) {
                    Field* current_field = fid->second;
                    // interpolate for every dimension
                    for (sid_t ida = 0; ida < current_field->lda(); ida++) {
                        // get the pointers
                        interp->Interpolate(dlvl, shift, &mem_src, block_out->data(current_field, ida), &mem_trg, block_in->data(current_field, ida));
                    }
                }
            }
        }
    }

    // deallocate the leaving blocks
    for (int id = 0; id < num_outgoing; id++) {
        qdrt_t*    quad  = outgoing[id];
        GridBlock* block = *(reinterpret_cast<GridBlock**>(quad->p.user_data));
        // delete the block, the fields are destroyed in the destructor
        delete (block);
    }

    // m_log("exit interpolate callback");
    m_profStop(grid->profiler(), "cback_interpolate");
    // m_log("exit prof callback");
    //-------------------------------------------------------------------------
    m_end;
}

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
    m_assert(grid->interp() != nullptr, "a Grid interpolator is needed");

    // get needed grid info
    auto                  f_start = grid->FieldBegin();
    auto                  f_end   = grid->FieldEnd();
    p8est_connectivity_t* connect = forest->connectivity;

    // allocate the incomming blocks
    for (int id = 0; id < num_incoming; id++) {
        qdrt_t* quad = incoming[id];
        // get block informations and create it
        real_t xyz[3];
        p8est_qcoord_to_vertex(connect, which_tree, quad->x, quad->y, quad->z, xyz);
        real_t     len   = m_quad_len(quad->level);
        GridBlock* block = new GridBlock(len, xyz, quad->level);
        // store the block
        *(reinterpret_cast<GridBlock**>(quad->p.user_data)) = block;
        // for every field, we allocate the memory
        for (auto fid = f_start; fid != f_end; fid++) {
            // allocate the new field
            block->AddField(fid->second);
        }

        // lock or unlock the block given the recursive behavior or not
        if (grid->recursive_adapt()) {
            block->unlock();
        } else {
            block->lock();
        }
    }

    // deallocate the leaving blocks
    for (int id = 0; id < num_outgoing; id++) {
        qdrt_t*    quad  = outgoing[id];
        GridBlock* block = *(reinterpret_cast<GridBlock**>(quad->p.user_data));
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
    Grid*     grid  = reinterpret_cast<Grid*>(forest->user_pointer);
    Field*    field = reinterpret_cast<Field*>(grid->cback_criterion_field());
    SetValue* expr  = reinterpret_cast<SetValue*>(grid->cback_interpolate_ptr());
    
    m_assert(grid->interp() != nullptr, "a Grid interpolator is needed");
    m_assert(expr->do_ghost(), "the SetValue object must set the ghost values");

    // get needed grid info
    auto                  f_start = grid->FieldBegin();
    auto                  f_end   = grid->FieldEnd();
    p8est_connectivity_t* connect = forest->connectivity;

    // allocate the incomming blocks
    for (int id = 0; id < num_incoming; id++) {
        qdrt_t* quad = incoming[id];
        // get block informations and create it
        real_t xyz[3];
        p8est_qcoord_to_vertex(connect, which_tree, quad->x, quad->y, quad->z, xyz);
        real_t     len   = m_quad_len(quad->level);
        GridBlock* block = new GridBlock(len, xyz, quad->level);
        // store the block
        *(reinterpret_cast<GridBlock**>(quad->p.user_data)) = block;
        // for every field, we allocate the memory
        for (auto fid = f_start; fid != f_end; fid++) {
            // allocate the new field
            block->AddField(fid->second);
        }
        // fill the new block with the analytical value
        expr->FillGridBlock(nullptr, block, field);

        // lock or unlock the block given the recursive behavior or not
        if (grid->recursive_adapt()) {
            block->unlock();
        } else {
            block->lock();
        }
    }

    // deallocate the leaving blocks
    for (int id = 0; id < num_outgoing; id++) {
        qdrt_t*    quad  = outgoing[id];
        GridBlock* block = *(reinterpret_cast<GridBlock**>(quad->p.user_data));
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
//     real_t     len      = m_quad_len(quad->level);
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