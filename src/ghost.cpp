#include "ghost.hpp"

#include "murphy.hpp"
#include "omp.h"
#include "wavelet.hpp"

// from paper p4est, table 1, S vector
static const sid_t edge2face[12][2]  = {{2, 4}, {3, 4}, {2, 5}, {3, 5}, {0, 4}, {1, 4}, {0, 5}, {1, 5}, {0, 2}, {1, 2}, {0, 3}, {1, 3}};
static const sid_t corner2face[8][3] = {{0, 2, 4}, {1, 2, 4}, {0, 3, 4}, {1, 3, 4}, {0, 2, 5}, {1, 2, 5}, {0, 3, 5}, {1, 3, 5}};

static const int facelimit[4] = {0, 24, 120, 144};
static const int edgelimit[4] = {0, 24, 72, 96};

Ghost::Ghost(Grid* grid,Interpolator* interpolator) {
    m_begin;
    m_assert(grid->is_mesh_valid(), "the mesh needs to be valid before entering here");
    //-------------------------------------------------------------------------
    grid_           = grid;
    interpolator_   = interpolator;
    p8est_t* forest = grid_->forest();
    // allocate the lists
    block_sibling_ = (list<GhostBlock*>**)m_calloc(forest->local_num_quadrants * sizeof(list<GhostBlock*>*));
    ghost_sibling_ = (list<GhostBlock*>**)m_calloc(forest->local_num_quadrants * sizeof(list<GhostBlock*>*));
    block_parent_  = (list<GhostBlock*>**)m_calloc(forest->local_num_quadrants * sizeof(list<GhostBlock*>*));
    ghost_parent_  = (list<GhostBlock*>**)m_calloc(forest->local_num_quadrants * sizeof(list<GhostBlock*>*));
    phys_          = (list<PhysBlock*>**)m_calloc(forest->local_num_quadrants * sizeof(list<PhysBlock*>*));

    // init the lists
    for (int ib = 0; ib < forest->local_num_quadrants; ib++) {
        // purge everything
        block_sibling_[ib] = new list<GhostBlock*>();
        ghost_sibling_[ib] = new list<GhostBlock*>();
        block_parent_[ib]  = new list<GhostBlock*>();
        ghost_parent_[ib]  = new list<GhostBlock*>();
        phys_[ib]          = new list<PhysBlock*>();
    }

    // call the simple operator to init the lists
    OperatorS::operator()(grid);

    // initialize the working coarse memory
    int nthreads = omp_get_max_threads();
    coarse_tmp_   = (real_p*)m_calloc(sizeof(real_p) * nthreads);
    for (int it = 0; it < nthreads; it++) {
        coarse_tmp_[it] =(real_t*) m_calloc(sizeof(real_t) * M_CLEN * M_CLEN * M_CLEN);
    }
    //-------------------------------------------------------------------------
    m_end;
}

Ghost::~Ghost() {
    m_begin;
    //-------------------------------------------------------------------------
    p8est_t* forest = grid_->forest();
    // clear the lists
    for (int ib = 0; ib < forest->local_num_quadrants; ib++) {
        // free the blocks
        for (auto biter = block_sibling_[ib]->begin(); biter != block_sibling_[ib]->end(); biter++) {
            delete (*biter);
        }
        for (auto biter = ghost_sibling_[ib]->begin(); biter != ghost_sibling_[ib]->end(); biter++) {
            delete (*biter);
        }
        for (auto biter = block_parent_[ib]->begin(); biter != block_parent_[ib]->end(); biter++) {
            delete (*biter);
        }
        for (auto biter = ghost_parent_[ib]->begin(); biter != ghost_parent_[ib]->end(); biter++) {
            delete (*biter);
        }
        for (auto piter = phys_[ib]->begin(); piter != phys_[ib]->end(); piter++) {
            delete (*piter);
        }
        // purge everything
        delete (block_sibling_[ib]);
        delete (ghost_sibling_[ib]);
        delete (block_parent_[ib]);
        delete (ghost_parent_[ib]);
        delete (phys_[ib]);
    }
    m_free(block_sibling_);
    m_free(ghost_sibling_);
    m_free(block_parent_);
    m_free(ghost_parent_);
    m_free(phys_);

    // free the temp memory
    int nthreads = omp_get_max_threads();
    for (int it = 0; it < nthreads; it++) {
        m_free(coarse_tmp_[it]);
    }
    m_free(coarse_tmp_);

    // free the memory
    m_free(mirrors_);
    m_free(ghosts_);
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief takes ghost values among neighbors and update the block's values
 * 
 * @param field 
 * @param current_ida 
 */
void Ghost::pull(Field* field){
    m_begin;
    //-------------------------------------------------------------------------
    
    for(sid_t ida = 0; ida<3; ida++){
        // setup the current ida
        ida_ = ida;
        // send the ghost


        // receive the ghosts
        
        // get a ghost copy
        OperatorF::operator()(grid_,field);
    }
    // set the ghost fields as ready
    field->ghost_status(true);
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief initialize the ghost data for a given block
 * 
 * @param qid 
 * @param block 
 */
void Ghost::apply(const qid_t* qid, GridBlock* block) {
    m_begin;
    //-------------------------------------------------------------------------
    // get the current lists
    list<GhostBlock*>* bsibling = block_sibling_[qid->cid];
    list<GhostBlock*>* gsibling = ghost_sibling_[qid->cid];
    list<GhostBlock*>* bparent  = block_parent_[qid->cid];
    list<GhostBlock*>* gparent  = ghost_parent_[qid->cid];
    list<PhysBlock*>*  phys     = phys_[qid->cid];

    // //----------------------------------
    // // temporary sc array used to get the ghosts
    sc_array_t*    ngh_quad = sc_array_new(sizeof(qdrt_t*));
    sc_array_t*    ngh_enc  = sc_array_new(sizeof(int));
    p8est_t*       forest   = grid_->forest();
    p8est_mesh_t*  mesh     = grid_->mesh();
    p8est_ghost_t* ghost    = grid_->ghost();

    for (sid_t ibidule = 0; ibidule < M_NNEIGHBOR; ibidule++) {
        // reset the quadrant and encoding stuff
        sc_array_reset(ngh_quad);
        sc_array_reset(ngh_enc);
        // sc_array_reset(ngh_qid);
        // get the neighboring quadrant
        p8est_mesh_get_neighbors(forest, ghost, mesh, qid->cid, ibidule, ngh_quad, ngh_enc, NULL);
        // decode the status and count the ghosts
        const size_t nghosts = ngh_enc->elem_count;
        //---------------------------------------------------------------------
        // we do the physics
        if (nghosts == 0) {
            sid_t isphys[3] = {0, 0, 0};
            // we only apply the physics to entire faces
            if (ibidule < 6) {
                phys->push_back(new PhysBlock(ibidule, block));
            }
            // else, the edges and corners will be filled through the face
        }
        //---------------------------------------------------------------------
        // this is a real block or a ghost
        for (int nid = 0; nid < nghosts; nid++) {
            const int  status  = *(ngh_enc->array + nid * sizeof(int));
            const bool isghost = (status < 0);
            qdrt_t*    nghq    = *((qdrt_t**)sc_array_index_int(ngh_quad, nid));

            // create a ghost block and store it given the gap in level
            if (!isghost) {
                GhostBlock* gb = new GhostBlock(block, nghq);
                if (gb->dlvl() >= 0) {
                    bsibling->push_back(gb);
                } else {
                    bparent->push_back(gb);
                }
            } else {
                // get the position of the tree, used to compute intersection between blocks
                p4est_topidx_t ngh_tree_id = nghq->p.piggy3.which_tree;
                real_t         ngh_tree_offset[3];
                p8est_qcoord_to_vertex(grid_->connect(), ngh_tree_id, 0, 0, 0, ngh_tree_offset);

                real_p      data = ghosts_ + qid->cid * M_NGHOST;
                GhostBlock* gb   = new GhostBlock(block, nghq, ngh_tree_offset, data);
                if (gb->dlvl() >= 0) {
                    gsibling->push_back(gb);
                } else {
                    gparent->push_back(gb);
                }
            }
        }
    }
    sc_array_destroy(ngh_quad);
    sc_array_destroy(ngh_enc);
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief computes the ghost points for a block given a field and the direction ida_
 * 
 * @param qid 
 * @param block 
 * @param fid 
 */
void Ghost::apply(const qid_t* qid, GridBlock* cur_block, Field* fid) {
    m_begin;
    //-------------------------------------------------------------------------
    // get the working direction given the thread
    const int ithread = omp_get_thread_num();
    real_p tmp = coarse_tmp_[ithread];

    // get a subblock describing the ghost memory
    lid_t     ghost_start[3] = {0, 0, 0};
    lid_t     ghost_end[3] = {M_N, M_N, M_N};
    SubBlock* ghost_subblock = new SubBlock(0, M_N, ghost_start, ghost_end);

    // determine if we have to use the coarse representationun
    const bool do_coarse = (block_parent_[qid->cid]->size() + ghost_parent_[qid->cid]->size()) > 0;

    for (auto biter = block_sibling_[qid->cid]->begin(); biter != block_sibling_[qid->cid]->end(); biter++) {
        GhostBlock* gblock = (*biter);
        // get the current blocks
        GridBlock* ngh_block = gblock->block_src();
        // get the memory pointers
        real_p data_src = ngh_block->data(fid,ida_);
        real_p data_trg = cur_block->data(fid,ida_);
        // launch the interpolation
        interpolator_->interpolate(gblock->dlvl(),gblock->shift(),ngh_block,data_src,gblock,data_trg);
    }
    for (auto biter = ghost_sibling_[qid->cid]->begin(); biter != ghost_sibling_[qid->cid]->end(); biter++) {
        GhostBlock* gblock = (*biter);
        // get the memory pointers for the ghost, where the send/recv call put the meaningfull info
        real_p data_src = gblock->data_src();
        real_p data_trg = cur_block->data(fid,ida_);
        // launch the interpolation
        interpolator_->interpolate(gblock->dlvl(),gblock->shift(),ghost_subblock,data_src,gblock,data_trg);
    }
    for (auto biter = block_parent_[qid->cid]->begin(); biter != block_parent_[qid->cid]->end(); biter++) {
        GhostBlock* gblock = (*biter);
    }
    for (auto biter = ghost_parent_[qid->cid]->begin(); biter != ghost_parent_[qid->cid]->end(); biter++) {
        GhostBlock* gblock = (*biter);
    }
    for (auto piter = phys_[qid->cid]->begin(); piter != phys_[qid->cid]->end(); piter++) {
        PhysBlock* gblock = (*piter);
    }


    delete (ghost_subblock);

    //-------------------------------------------------------------------------
    m_end;
}