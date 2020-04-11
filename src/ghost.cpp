#include "ghost.hpp"

#include "murphy.hpp"

// from paper p4est, table 1, S vector
static const sid_t edge2face[12][2]  = {{2, 4}, {3, 4}, {2, 5}, {3, 5}, {0, 4}, {1, 4}, {0, 5}, {1, 5}, {0, 2}, {1, 2}, {0, 3}, {1, 3}};
static const sid_t corner2face[8][3] = {{0, 2, 4}, {1, 2, 4}, {0, 3, 4}, {1, 3, 4}, {0, 2, 5}, {1, 2, 5}, {0, 3, 5}, {1, 3, 5}};

static const int facelimit[4] = {0, 24, 120, 144};
static const int edgelimit[4] = {0, 24, 72, 96};

Ghost::Ghost(Grid* grid) {
    m_begin;
    m_assert(grid->is_mesh_valid(),"the mesh needs to be valid before entering here");
    //-------------------------------------------------------------------------
    grid_           = grid;
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
        for (list<GhostBlock*>::iterator biter = block_sibling_[ib]->begin(); biter != block_sibling_[ib]->end(); biter++) {
            delete (*biter);
        }
        for (list<GhostBlock*>::iterator biter = ghost_sibling_[ib]->begin(); biter != ghost_sibling_[ib]->end(); biter++) {
            delete (*biter);
        }
        for (list<GhostBlock*>::iterator biter = block_parent_[ib]->begin(); biter != block_parent_[ib]->end(); biter++) {
            delete (*biter);
        }
        for (list<GhostBlock*>::iterator biter = ghost_parent_[ib]->begin(); biter != ghost_parent_[ib]->end(); biter++) {
            delete (*biter);
        }
        for (list<PhysBlock*>::iterator  piter = phys_[ib]->begin(); piter != phys_[ib]->end(); piter++) {
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

    // free the memory
    m_free(mirrors_);
    m_free(ghosts_);
    //-------------------------------------------------------------------------
    m_end;
}

void Ghost::apply(const qid_t* qid, GridBlock* block)  {
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
    // sc_array_t*    ngh_qid  = sc_array_new(sizeof(int));
    p8est_t*       forest   = grid_->forest();
    p8est_mesh_t*  mesh     = grid_->mesh();
    p8est_ghost_t* ghost    = grid_->ghost();

    for (int ibidule = 0; ibidule < M_NNEIGHBOR; ibidule++) {
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
    // sc_array_destroy(ngh_qid);

    //-------------------------------------------------------------------------
    m_end;
}



void Ghost::apply(const qid_t* qid, GridBlock* block, const Field* fid){
    m_begin;
    //-------------------------------------------------------------------------
    printf("kikouu");
    //-------------------------------------------------------------------------
    m_end;
}