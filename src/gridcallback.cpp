
#include "gridcallback.hpp"

#include <p8est_bits.h>

#include "gridblock.hpp"
#include "grid.hpp"

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

    real_t len        = m_quad_len(quad->level);
    // the user data points to the data defined by the grid = GridBlock*
    *(reinterpret_cast<GridBlock**>(quad->p.user_data)) = new GridBlock(len, xyz, quad->level);
    //-------------------------------------------------------------------------
    m_end;
}

void cback_DestroyBlock(p8est_iter_volume_info_t* info, void* user_data) {
    m_begin;
    //-------------------------------------------------------------------------
    p8est_quadrant_t* quad = info->quad;
    GridBlock* block = *(reinterpret_cast<GridBlock**>(quad->p.user_data)) ;
    delete (block);
    *(reinterpret_cast<GridBlock**>(quad->p.user_data)) = nullptr;
    //-------------------------------------------------------------------------
    m_end;
}



int cback_Yes(p8est_t *forest, p4est_topidx_t which_tree, qdrt_t *quadrant) {
    m_begin;
    //-------------------------------------------------------------------------
    return (true);
    //-------------------------------------------------------------------------
    m_end;
}
int cback_Yes(p8est_t *forest, p4est_topidx_t which_tree, qdrt_t *quadrant[]) {
    m_begin;
    //-------------------------------------------------------------------------
    return (true);
    //-------------------------------------------------------------------------
    m_end;
}

// this is refinement
int cback_Wavelet(p8est_t *forest, p4est_topidx_t which_tree, qdrt_t *quadrant) {
    m_begin;
    //-------------------------------------------------------------------------
    Grid*         grid   = reinterpret_cast<Grid*>(forest->user_pointer);
    GridBlock*    block   = *(reinterpret_cast<GridBlock**>(quadrant->p.user_data));
    Interpolator* interp = grid->interp();
    // get the field and check each dimension
    bool   refine = false;
    Field* fid    = grid->tmp_field();
    for (int ida = 0; ida < fid->lda(); ida++) {
        m_verb("analysing the details on block for field %s, dim = %d",fid->name().c_str(),ida);
        real_p data = block->data(fid, ida);
        real_t norm = interp->Criterion(block, data);
        // refine if the norm is bigger
        refine = (norm > grid->rtol());
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

// this is coarsening
int cback_Wavelet(p8est_t *forest, p4est_topidx_t which_tree, qdrt_t *quadrant[]) {
    m_begin;
    //-------------------------------------------------------------------------
    Grid* grid = reinterpret_cast<Grid*>(forest->user_pointer);
    Interpolator* interp = grid->interp();
    // for each of the children
    bool   coarsen = false;
    Field* fid     = grid->tmp_field();
    for (int id = 0; id < P8EST_CHILDREN; id++) {
        GridBlock* block = *(reinterpret_cast<GridBlock**>(quadrant[id]->p.user_data));
        for (int ida = 0; ida < fid->lda(); ida++) {
            real_p data = block->data(fid, ida);
            real_t norm = interp->Criterion(block, data);
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

void cback_Interpolate(p8est_t* forest, p4est_topidx_t which_tree, int num_outgoing, qdrt_t* outgoing[], int num_incoming, qdrt_t* incoming[]) {
    m_begin;
    m_assert(num_incoming == 1 || num_outgoing == 1,"we have either to compress or to refine");
    m_assert(num_incoming == P8EST_CHILDREN || num_outgoing == P8EST_CHILDREN,"the number of replacing blocks has to be the number of children");
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

    // get the working field
    // Field* field = grid->working_callback_field();
    
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
    const sid_t dlvl = (num_outgoing==1) ? -1 : 1;

    // create an empty SubBlock representing the valid GP
    lid_t     full_start[3]  = {0, 0, 0};
    lid_t     full_end[3]    = {0, 0, 0};
    SubBlock* mem_block = new SubBlock(M_GS,M_STRIDE,full_start,full_end);

    // do the required iterpolation
    for(sid_t iout=0; iout<num_outgoing; iout++){
        GridBlock* block_out = *(reinterpret_cast<GridBlock**>(outgoing[iout]->p.user_data));

        for (sid_t iin = 0; iin < num_incoming; iin++) {
            GridBlock* block_in = *(reinterpret_cast<GridBlock**>(incoming[iout]->p.user_data));
            
            
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
                        interp->Interpolate(dlvl, shift, mem_block, block_out->data(current_field, ida), block_in, block_in->data(current_field, ida));
                    }
                }
            }
            // if we coarsen
            else if (num_incoming == 1) {
                // get the shift
                lid_t shift[3];
                shift[0] = - M_N * ((childid % 2));      // x corner index = (ic%4)%2
                shift[1] = - M_N * ((childid % 4) / 2);  // y corner index
                shift[2] = - M_N * ((childid / 4));      // z corner index
                // we create a subblock with the correct memory representation for the source
                // no nned to redefine the source zone, we can use the standard 0 -> M_N
                full_start[0] = M_HN * ((childid % 2));      // x corner index = (ic%4)%2
                full_start[1] = M_HN * ((childid % 4) / 2);  // y corner index
                full_start[2] = M_HN * ((childid / 4));      // z corner index
                for (int id = 0; id < 3; id++) {
                    full_end[id]   = full_start[id] + M_HN;
                }
                mem_block->Reset(M_GS, M_STRIDE, full_start, full_end);

                // for every field, we interpolate it
                for (auto fid = f_start; fid != f_end; fid++) {
                    Field* current_field = fid->second;
                    // interpolate for every dimension
                    for (sid_t ida = 0; ida < current_field->lda(); ida++) {
                        // get the pointers
                        interp->Interpolate(dlvl, shift, block_out, block_out->data(current_field, ida), mem_block, block_in->data(current_field, ida));
                    }
                }
            }
        }
    }

    delete(mem_block);

    // deallocate the leaving blocks
    for (int id = 0; id < num_outgoing; id++) {
        qdrt_t*    quad  = outgoing[id];
        GridBlock* block = *(reinterpret_cast<GridBlock**>(quad->p.user_data));
        // delete the block, the fields are destroyed in the destructor
        delete(block);
    }
    //-------------------------------------------------------------------------
    m_end;
}