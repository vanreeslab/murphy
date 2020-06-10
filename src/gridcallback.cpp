#include "gridcallback.hpp"

#include <p8est_bits.h>

#include <list>

#include "grid.hpp"
#include "gridblock.hpp"
#include "patch.hpp"
#include "mgfamily.hpp"
#include "multigrid.hpp"

/**
 * @brief initiate a new block and store its adress in the p4est quad
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
    list<Patch>* patches = reinterpret_cast<list<Patch>*>(grid->tmp_ptr());
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
    list<Patch>* patches = reinterpret_cast<list<Patch>*>(grid->tmp_ptr());

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
 * @brief reply that we should refine the block if the output of @ref Interpolator::Criterion() is bigger than the tolerance (we use the interpolator from the @ref Grid)
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
    GridBlock*    block  = *(reinterpret_cast<GridBlock**>(quadrant->p.user_data));
    Interpolator* interp = grid->interp();
    // get the field and check each dimension
    bool   refine = false;
    Field* fid    = reinterpret_cast<Field*>(grid->tmp_ptr());
    for (int ida = 0; ida < fid->lda(); ida++) {
        real_p data = block->data(fid, ida);
        real_t norm = interp->Criterion(block, data,grid->mem_pool());
        // refine if the norm is bigger
        refine = (norm > grid->rtol());
        m_verb("refine? %e vs %e", norm, grid->rtol());
        // refine the whole grid if one dimension needs to be refined
        if (refine) {
            break;
        }
    }
    // refine if the criterion is bigger than the tolerance
    return refine;
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief reply that we should coarsen the block if the output of the @ref Interpolator::Criterion() is bigger than the tolerance (we use the interpolator from the @ref Grid)
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
    Interpolator* interp = grid->interp();
    // for each of the children
    bool   coarsen = false;
    Field* fid    = reinterpret_cast<Field*>(grid->tmp_ptr());
    for (int id = 0; id < P8EST_CHILDREN; id++) {
        GridBlock* block = *(reinterpret_cast<GridBlock**>(quadrant[id]->p.user_data));
        for (int ida = 0; ida < fid->lda(); ida++) {
            real_p data = block->data(fid, ida);
            real_t norm = interp->Criterion(block, data,grid->mem_pool());
            // coarsen if the norm is bigger
            coarsen = (norm < grid->ctol());
            // coarsen the whole grid if one dimension needs to be coarsened
            if (coarsen) {
                break;
            }
        }
        // coarsen the whole grid if one block needs to be coarsened
        if (coarsen) {
            break;
        }
    }
    // refine if the criterion is bigger than the tolerance
    return coarsen;
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief Interpolate one block to his children or 8 children to their parent, using the Interpolator of the @ref Grid.
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

    // get needed grid info
    auto                  f_start = grid->FieldBegin();
    auto                  f_end   = grid->FieldEnd();
    Interpolator*         interp  = grid->interp();
    p8est_connectivity_t* connect = forest->connectivity;

    if (grid->HasProfiler()) {
        grid->profiler()->Start("cback_interpolate");
    }

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

    // create an empty SubBlock representing the valid GP
    lid_t     full_start[3] = {0, 0, 0};
    lid_t     full_end[3]   = {0, 0, 0};
    SubBlock* mem_block     = new SubBlock(M_GS, M_STRIDE, full_start, full_end);

    // do the required iterpolation
    for (sid_t iout = 0; iout < num_outgoing; iout++) {
        GridBlock* block_out = *(reinterpret_cast<GridBlock**>(outgoing[iout]->p.user_data));

        for (sid_t iin = 0; iin < num_incoming; iin++) {
            GridBlock* block_in = *(reinterpret_cast<GridBlock**>(incoming[iin]->p.user_data));

            // get the correct child id: if 1 is going out, the children are the incoming
            int childid = p8est_quadrant_child_id((num_outgoing == 1) ? incoming[iin] : outgoing[iout]);

            // if we refine
            if (num_outgoing == 1) {
                // get the shift given the child id
                lid_t shift[3];
                shift[0] = M_HN * ((childid % 2));      // x corner index = (ic%4)%2
                shift[1] = M_HN * ((childid % 4) / 2);  // y corner index
                shift[2] = M_HN * ((childid / 4));      // z corner index
                // we create a subblock with the correct memory representation for the target
                // no nned to redefine the target zone, we can use the standard 0 -> M_N
                for (int id = 0; id < 3; id++) {
                    full_start[id] = shift[id] - M_GS;
                    full_end[id]   = shift[id] + M_HN + M_GS;
                }
                mem_block->Reset(M_GS, M_STRIDE, full_start, full_end);

                // for every field, we interpolate it
                for (auto fid = f_start; fid != f_end; fid++) {
                    Field* current_field = fid->second;
                    // interpolate for every dimension
                    for (sid_t ida = 0; ida < current_field->lda(); ida++) {
                        // get the pointers
                        interp->Interpolate(dlvl, shift, mem_block, block_out->data(current_field, ida), block_in, block_in->data(current_field, ida), grid->mem_pool());
                    }
                }
            } else if (num_incoming == 1) {
                // get the shift
                lid_t shift[3];
                shift[0] = -M_N * ((childid % 2));      // x corner index = (ic%4)%2
                shift[1] = -M_N * ((childid % 4) / 2);  // y corner index
                shift[2] = -M_N * ((childid / 4));      // z corner index
                // we create a subblock with the correct memory representation for the source
                // no nned to redefine the source zone, we can use the standard 0 -> M_N
                full_start[0] = M_HN * ((childid % 2));      // x corner index = (ic%4)%2
                full_start[1] = M_HN * ((childid % 4) / 2);  // y corner index
                full_start[2] = M_HN * ((childid / 4));      // z corner index
                for (int id = 0; id < 3; id++) {
                    full_end[id] = full_start[id] + M_HN;
                }
                mem_block->Reset(M_GS, M_STRIDE, full_start, full_end);

                // for every field, we interpolate it
                for (auto fid = f_start; fid != f_end; fid++) {
                    Field* current_field = fid->second;
                    // interpolate for every dimension
                    for (sid_t ida = 0; ida < current_field->lda(); ida++) {
                        // get the pointers
                        interp->Interpolate(dlvl, shift, block_out, block_out->data(current_field, ida), mem_block, block_in->data(current_field, ida), grid->mem_pool());
                    }
                }
            }
        }
    }

    delete (mem_block);

    // deallocate the leaving blocks
    for (int id = 0; id < num_outgoing; id++) {
        qdrt_t*    quad  = outgoing[id];
        GridBlock* block = *(reinterpret_cast<GridBlock**>(quad->p.user_data));
        // delete the block, the fields are destroyed in the destructor
        delete (block);
    }
    if (grid->HasProfiler()) {
        grid->profiler()->Stop("cback_interpolate");
    }
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

void cback_MGCreateFamilly(p8est_t* forest, p4est_topidx_t which_tree, int num_outgoing, qdrt_t* outgoing[], int num_incoming, qdrt_t* incoming[]) {
    m_begin;
    m_assert(num_outgoing == P8EST_CHILDREN, "this function is only called when doing the refinement");
    m_assert(num_incoming == 1, "this function is only called when doing the refinement");
    //-------------------------------------------------------------------------
    p8est_connectivity_t* connect = forest->connectivity;
    // retrieve the multigrid
    Multigrid* grid = reinterpret_cast<Multigrid*>(forest->user_pointer);

    // get the new block informations and create it
    qdrt_t* quad = incoming[0];
    real_t xyz[3];
    p8est_qcoord_to_vertex(connect, which_tree, quad->x, quad->y, quad->z, xyz);
    real_t     len      = m_quad_len(quad->level);
    GridBlock* parent = new GridBlock(len, xyz, quad->level);
    // allocate the correct fields
    parent->AddFields(grid->map_fields());
    // store the new block
    *(reinterpret_cast<GridBlock**>(quad->p.user_data)) = parent;

    // get the leaving blocks
    GridBlock* children[P8EST_CHILDREN];
    for(sid_t ic=0; ic<P8EST_CHILDREN; ic++){
        children[ic] = *(reinterpret_cast<GridBlock**>(outgoing[ic]->p.user_data));
    }
    
    // bind the family together
    MGFamily* family = grid->curr_family();
    family->AddMembers(parent,children);
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief always reply yes if p4est ask if we should coarsen the block
 */
int cback_Level(p8est_t* forest, p4est_topidx_t which_tree, qdrt_t* quadrant[]) {
    m_begin;
    //-------------------------------------------------------------------------
    Multigrid* grid = reinterpret_cast<Multigrid*>(forest->user_pointer);
    sid_t target_level = grid->curr_level();
    return (quadrant[0]->level > target_level);
    //-------------------------------------------------------------------------
    m_end;
}