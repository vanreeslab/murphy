#include "ghost.hpp"

#include "boundary.hpp"
#include "murphy.hpp"
#include "omp.h"
#include "wavelet.hpp"

// // from paper p4est, table 1, S vector
// static const sid_t edge2face[12][2]  = {{2, 4}, {3, 4}, {2, 5}, {3, 5}, {0, 4}, {1, 4}, {0, 5}, {1, 5}, {0, 2}, {1, 2}, {0, 3}, {1, 3}};
// static const sid_t corner2face[8][3] = {{0, 2, 4}, {1, 2, 4}, {0, 3, 4}, {1, 3, 4}, {0, 2, 5}, {1, 2, 5}, {0, 3, 5}, {1, 3, 5}};

// static const int facelimit[4] = {0, 24, 120, 144};
// static const int edgelimit[4] = {0, 24, 72, 96};

Ghost::Ghost(ForestGrid* grid) {
    m_begin;
    m_assert(grid->is_mesh_valid(), "the mesh needs to be valid before entering here");
    //-------------------------------------------------------------------------
    // get the important pointers
    grid_              = grid;
    p8est_mesh_t* mesh = grid->mesh();

    // allocate the lists
    block_sibling_ = (list<GhostBlock*>**)m_calloc(mesh->local_num_quadrants * sizeof(list<GhostBlock*>*));
    ghost_sibling_ = (list<GhostBlock*>**)m_calloc(mesh->local_num_quadrants * sizeof(list<GhostBlock*>*));
    block_parent_  = (list<GhostBlock*>**)m_calloc(mesh->local_num_quadrants * sizeof(list<GhostBlock*>*));
    ghost_parent_  = (list<GhostBlock*>**)m_calloc(mesh->local_num_quadrants * sizeof(list<GhostBlock*>*));
    phys_          = (list<PhysBlock*>**)m_calloc(mesh->local_num_quadrants * sizeof(list<PhysBlock*>*));

    // init the lists
    for (int ib = 0; ib < mesh->local_num_quadrants; ib++) {
        // purge everything
        block_sibling_[ib] = new list<GhostBlock*>();
        ghost_sibling_[ib] = new list<GhostBlock*>();
        block_parent_[ib]  = new list<GhostBlock*>();
        ghost_parent_[ib]  = new list<GhostBlock*>();
        phys_[ib]          = new list<PhysBlock*>();
    }

    // initialize the communications and the mirrors, ghosts arrays
    InitComm_();

    // call the simple operator to init the lists, reset the counter to 0 (needed to count the ghosts)
    OperatorS::operator()(grid_);

    // initialize the working coarse memory
    int nthreads = omp_get_max_threads();
    coarse_tmp_  = (real_p*)m_calloc(sizeof(real_p) * nthreads);
    for (int it = 0; it < nthreads; it++) {
        coarse_tmp_[it] = (real_t*)m_calloc(sizeof(real_t) * M_CLEN * M_CLEN * M_CLEN);
    }
    //-------------------------------------------------------------------------
    m_log("ghost ready to exchange");
    m_end;
}

Ghost::~Ghost() {
    m_begin;
    m_assert(grid_->is_mesh_valid(), "the mesh needs to be valid before entering here");
    //-------------------------------------------------------------------------
    p8est_mesh_t* mesh = grid_->mesh();
    // clear the lists
    for (lid_t ib = 0; ib < mesh->local_num_quadrants; ib++) {
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

    // free the requests
    for (int ir = 0; ir < n_send_request_; ir++) {
        MPI_Request_free(mirror_send_ + ir);
    }
    m_free(mirror_send_);
    for (int ir = 0; ir < n_recv_request_; ir++) {
        MPI_Request_free(ghost_recv_ + ir);
    }
    m_free(ghost_recv_);

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

void Ghost::PushToMirror(Field* field, sid_t ida) {
    m_begin;
    m_assert(grid_->is_mesh_valid(), "the mesh needs to be valid before entering here");
    //-------------------------------------------------------------------------
    // store the current interpolator and dimension
    ida_ = ida;
    // loop on the mirrors and copy the values
    LoopOnMirrorBlock_(&Ghost::PushToMirror_, field);
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief starts the send requests of the mirror
 * 
 * @warning we do not start the reception requests because we are not sure the buffers are available
 * 
 */
void Ghost::MirrorToGhostSend() {
    m_begin;
    //-------------------------------------------------------------------------
    if (n_send_request_ > 0) {
        MPI_Startall(n_send_request_, mirror_send_);
    }
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief starts the reception requests, end the sending request and end the reception requests
 * 
 */
void Ghost::MirrorToGhostRecv() {
    m_begin;
    //-------------------------------------------------------------------------
    if (n_recv_request_ > 0) {
        MPI_Startall(n_recv_request_, ghost_recv_);
    }
    if (n_send_request_ > 0) {
        MPI_Waitall(n_send_request_, mirror_send_, MPI_STATUSES_IGNORE);
    }
    if (n_recv_request_ > 0) {
        MPI_Waitall(n_recv_request_, ghost_recv_, MPI_STATUSES_IGNORE);
    }
    //-------------------------------------------------------------------------
    m_end;
}

void Ghost::PullFromGhost(Field* field, sid_t ida, Interpolator* interp) {
    m_begin;
    m_assert(grid_->is_mesh_valid(), "the mesh needs to be valid before entering here");
    //-------------------------------------------------------------------------
    // store the current interpolator and dimension
    interp_ = interp;
    ida_    = ida;
    // interpolate
    OperatorF::operator()(grid_, field);
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief initialize the ghost data for a given block
 * 
 * @param qid 
 * @param block 
 */
void Ghost::ApplyOperatorS(const qid_t* qid, GridBlock* block) {
    //-------------------------------------------------------------------------
    InitList_(qid, block);
    //-------------------------------------------------------------------------
}

/**
 * @brief computes the ghost points for a block, given a field and the direction ida_
 * 
 * @param qid 
 * @param block 
 * @param fid 
 */
void Ghost::ApplyOperatorF(const qid_t* qid, GridBlock* cur_block, Field* fid) {
    //-------------------------------------------------------------------------
    PullFromGhost_(qid, cur_block, fid);
    //-------------------------------------------------------------------------
}

void Ghost::InitComm_() {
    m_begin;
    //-------------------------------------------------------------------------
    int            mpi_size = grid_->mpisize();
    MPI_Comm       mpi_comm = grid_->mpicomm();
    p8est_ghost_t* ghost    = grid_->ghost();

    // this is the sum of mirrors to send, as a given mirror can be send to multiple ranks
    // this will not happen for a ghost
    n_mirror_to_send_ = 0;
    n_send_request_   = 0;
    n_recv_request_   = 0;
    for (int ir = 0; ir < mpi_size; ir++) {
        // send
        lid_t send_first  = ghost->mirror_proc_offsets[ir];
        lid_t send_last   = ghost->mirror_proc_offsets[ir + 1];
        n_send_request_   = n_send_request_ + ((send_last - send_first) > 0 ? 1 : 0);
        n_mirror_to_send_ = n_mirror_to_send_ + (send_last - send_first);
        // recv
        lid_t recv_first = ghost->proc_offsets[ir];
        lid_t recv_last  = ghost->proc_offsets[ir + 1];
        n_recv_request_ += (recv_last - recv_first) > 0 ? 1 : 0;
    }
    // allocate the request arrays
    mirror_send_ = (MPI_Request*)m_calloc(n_send_request_ * sizeof(MPI_Request));
    ghost_recv_  = (MPI_Request*)m_calloc(n_recv_request_ * sizeof(MPI_Request));

    // allocate the mirror and ghost arrayss
    mirrors_to_local_ = (lid_t*)m_calloc(sizeof(lid_t) * n_mirror_to_send_);
    mirrors_          = (real_t*)m_calloc(sizeof(real_t) * M_NGHOST * n_mirror_to_send_);
    ghosts_           = (real_t*)m_calloc(sizeof(real_t) * M_NGHOST * ghost->ghosts.elem_count);

    // allocate the requests
    lid_t  send_request_offset = 0;
    lid_t  recv_request_offset = 0;
    size_t mirror_offset       = 0;
    for (int ir = 0; ir < mpi_size; ir++) {
        // send
        // mirror proc offset tells how many to send
        // mirror proc mirror tells which mirror is to send
        lid_t send_first = ghost->mirror_proc_offsets[ir];
        lid_t send_last  = ghost->mirror_proc_offsets[ir + 1];
        lid_t block2send = send_last - send_first;
        if (block2send > 0) {
            // get the starting buffer adress and the size
            real_p send_buf  = mirrors_ + mirror_offset;
            size_t send_size = block2send * M_NGHOST;
            // init the send
            MPI_Send_init(send_buf, send_size, M_MPI_REAL, ir, P4EST_COMM_GHOST_EXCHANGE, mpi_comm, mirror_send_ + send_request_offset);
            // store the mirror id's
            for (lid_t ib = send_first; ib < send_last; ib++) {
                mirrors_to_local_[ib] = ghost->mirror_proc_mirrors[ib];
            }
            // upate the counters
            send_request_offset = send_request_offset + 1;
            mirror_offset       = mirror_offset + send_size;
        }

        // recv
        lid_t recv_first = ghost->proc_offsets[ir];
        lid_t recv_last  = ghost->proc_offsets[ir + 1];
        lid_t block2recv = recv_last - recv_first;
        if (block2recv > 0) {
            // get the starting buffer adress and the size
            real_p recv_buf  = ghosts_ + recv_first * M_NGHOST;
            int    recv_size = block2recv * M_NGHOST;
            // init the send
            MPI_Recv_init(recv_buf, recv_size, M_MPI_REAL, ir, P4EST_COMM_GHOST_EXCHANGE, mpi_comm, ghost_recv_ + recv_request_offset);
            // upate the counters
            recv_request_offset = recv_request_offset + 1;
        }
    }

    m_assert(send_request_offset == n_send_request_, "the two numbers have to match");
    m_assert(recv_request_offset == n_recv_request_, "the two numbers have to match");

    //-------------------------------------------------------------------------
    m_end;
}

void Ghost::InitList_(const qid_t* qid, GridBlock* block) {
    //-------------------------------------------------------------------------
    // get the current lists
    list<GhostBlock*>* bsibling = block_sibling_[qid->cid];
    list<GhostBlock*>* gsibling = ghost_sibling_[qid->cid];
    list<GhostBlock*>* bparent  = block_parent_[qid->cid];
    list<GhostBlock*>* gparent  = ghost_parent_[qid->cid];
    list<PhysBlock*>*  phys     = phys_[qid->cid];

    // //----------------------------------
    // // temporary sc array used to get the ghosts
    sc_array_t* ngh_quad;  // points to the quad
    sc_array_t* ngh_qid;   // give the ID of the quad or ghost
    sc_array_t* ngh_enc;   // get the status

#pragma omp critical
    {
        ngh_quad = sc_array_new(sizeof(qdrt_t*));
        ngh_enc  = sc_array_new(sizeof(int));
        ngh_qid  = sc_array_new(sizeof(int));
    }

    p8est_t*       forest = grid_->forest();
    p8est_mesh_t*  mesh   = grid_->mesh();
    p8est_ghost_t* ghost  = grid_->ghost();

    for (sid_t ibidule = 0; ibidule < M_NNEIGHBOR; ibidule++) {
#pragma omp critical
        {
            // get the neighboring quadrant
            sc_array_reset(ngh_quad);
            sc_array_reset(ngh_enc);
            sc_array_reset(ngh_qid);
            p8est_mesh_get_neighbors(grid_->forest(), grid_->ghost(), grid_->mesh(), qid->cid, ibidule, ngh_quad, ngh_enc, ngh_qid);
        }
        // decode the status and count the ghosts
        const size_t nghosts = ngh_enc->elem_count;
        //---------------------------------------------------------------------
        // we do the physics
        if (nghosts == 0) {
            sid_t isphys[3] = {0, 0, 0};
            // we only apply the physics to entire faces
            if (ibidule < 6) {
                PhysBlock* pb = new PhysBlock(ibidule, block);
#pragma omp critical
                phys->push_back(pb);
            }
            // else, the edges and corners will be filled through the face
        }
        //---------------------------------------------------------------------
        // this is a real block or a ghost
        for (int nid = 0; nid < nghosts; nid++) {
            const int  status  = *(ngh_enc->array + nid * sizeof(int));
            const bool isghost = (status < 0);
            qdrt_t*    nghq    = *((qdrt_t**)sc_array_index_int(ngh_quad, nid));

            // get the sign, i.e. the normal to the face, the edge of the corner we consider
            real_t sign[3];
            GhostGetSign(ibidule, sign);

            // get the position of the neighbor, as seen by me!!! may be different than the actual one if there is a periodic bc
            real_t ngh_pos[3];
            if (!isghost) {
                // cannot use the p8est function because the which_tree is not assigned, so we retrieve the position through the block
                GridBlock* ngh_block = reinterpret_cast<GridBlock*>(nghq->p.user_data);
                ngh_pos[0]           = ngh_block->xyz(0);
                ngh_pos[1]           = ngh_block->xyz(1);
                ngh_pos[2]           = ngh_block->xyz(2);
            } else {
                p8est_qcoord_to_vertex(grid_->connect(), nghq->p.piggy3.which_tree, nghq->x, nghq->y, nghq->z, ngh_pos);
            }
            for (sid_t id = 0; id < 3; id++) {
                ngh_pos[id] = block->xyz(id) + fmod(ngh_pos[id] - block->xyz(id) + grid_->domain_periodic(id)*sign[id] * grid_->domain_length(id), grid_->domain_length(id));
            }

            // create the new block
            GhostBlock* gb = new GhostBlock(block, nghq->level, ngh_pos);

            // associate the correct memory and push back
            if (!isghost) {
                GridBlock* ngh_block = reinterpret_cast<GridBlock*>(nghq->p.user_data);
                gb->block_src(ngh_block);

                if (gb->dlvl() >= 0) {
#pragma omp critical
                    bsibling->push_back(gb);
                } else {
#pragma omp critical
                    bparent->push_back(gb);
                }
                // m_verb("dbg: tree %d, quad %d; block detected for ibidule %d",qid->tid,qid->qid,ibidule);
            } else {
                lid_t  ighost = *(ngh_qid->array + nid * sizeof(int));
                real_p data   = ghosts_ + ighost * M_NGHOST;
                gb->data_src(data);

                if (gb->dlvl() >= 0) {
#pragma omp critical
                    gsibling->push_back(gb);
                } else {
#pragma omp critical
                    gparent->push_back(gb);
                }
                // m_verb("dbg: tree %d, quad %d; ghost detected for ibidule %d",qid->tid,qid->qid,ibidule);
            }
        }
    }

#pragma omp critical
    {
        sc_array_destroy(ngh_quad);
        sc_array_destroy(ngh_enc);
        sc_array_destroy(ngh_qid);
    }
    //-------------------------------------------------------------------------
}

void Ghost::PushToMirror_(const qid_t* qid, GridBlock* block, Field* fid) {
    m_assert(ida_ >= 0, "the current working dimension has to be correct");
    m_assert(ida_ < fid->lda(), "the current working dimension has to be correct");
    //-------------------------------------------------------------------------
    real_p mirror = mirrors_ + qid->cid * M_NGHOST;
    real_p data   = block->data(fid, ida_);

    for (int i2 = 0; i2 < M_N; i2++) {
        for (int i1 = 0; i1 < M_N; i1++) {
            for (int i0 = 0; i0 < M_N; i0++) {
                mirror[m_sidx(i0, i1, i2, 0, M_N)] = data[m_idx(i0, i1, i2)];
            }
        }
    }
    //-------------------------------------------------------------------------
}

/**
 * @brief given an updated @ref ghost_ arrays,
 * compute the ghosts points for the current block (field @ref fid, dimension @ref ida_)
 * 
 * @warning this function runs inside an OpenMP parallel region
 * 
 * @param qid the quadrant id, see @ref qid_t
 * @param cur_block the current grid block, see @ref GridBlock
 * @param fid the current field, see @ref Field
 */
void Ghost::PullFromGhost_(const qid_t* qid, GridBlock* cur_block, Field* fid) {
    //-------------------------------------------------------------------------
    // get the working direction given the thread
    const int ithread = omp_get_thread_num();
    real_p    tmp     = coarse_tmp_[ithread];

    // get a subblock describing the ghost memory
    lid_t     ghost_start[3] = {0, 0, 0};
    lid_t     ghost_end[3]   = {M_N, M_N, M_N};
    SubBlock* ghost_subblock = new SubBlock(0, M_N, ghost_start, ghost_end);

    lid_t     coarse_start[3] = {0, 0, 0};
    lid_t     coarse_end[3]   = {M_HN, M_HN, M_HN};
    SubBlock* coarse_subblock = new SubBlock(M_GS, M_CLEN, coarse_start, coarse_end);

    // determine if we have to use the coarse representationun
    const bool do_coarse = (block_parent_[qid->cid]->size() + ghost_parent_[qid->cid]->size()) > 0;
    // if so, reset the coarse info
    if (do_coarse) {
        memset(tmp, 0, M_CLEN * M_CLEN * M_CLEN * sizeof(real_t));
    }

    m_verb("dbg: tree %d, quad %d: %ld blocks, %ld ghosts, %ld phys", qid->tid, qid->qid, block_sibling_[qid->cid]->size(), ghost_sibling_[qid->cid]->size(), phys_[qid->cid]->size());

    //-------------------------------------------------------------------------
    // do the blocks, on my level or finer
    for (auto biter = block_sibling_[qid->cid]->begin(); biter != block_sibling_[qid->cid]->end(); biter++) {
        GhostBlock* gblock    = (*biter);
        GridBlock*  ngh_block = gblock->block_src();
        // memory details
        MemLayout* block_src = ngh_block;
        real_p     data_src  = ngh_block->data(fid, ida_);
        MemLayout* block_trg = gblock;
        real_p     data_trg  = cur_block->data(fid, ida_);
        // interpolate
        interp_->Interpolate(gblock->dlvl(), gblock->shift(), block_src, data_src, block_trg, data_trg);

        // we need to interpolate on the coarse version of myself as well
        if (do_coarse) {
            // set the coarse block to the correct position
            for (int id = 0; id < 3; id++) {
                coarse_start[id] = CoarseFromBlock(gblock->start(id));
                coarse_end[id]   = CoarseFromBlock(gblock->end(id));
            }
            coarse_subblock->Reset(M_GS, M_CLEN, coarse_start, coarse_end);
            // memory details
            MemLayout* block_src = ngh_block;
            real_p     data_src  = ngh_block->data(fid, ida_);
            MemLayout* block_trg = coarse_subblock;
            real_p     data_trg  = tmp + m_zeroidx(0, coarse_subblock);
            // interpolate, the level is 1 coarser and the shift is unchanged
            interp_->Interpolate(gblock->dlvl() + 1, gblock->shift(), block_src, data_src, block_trg, data_trg);
        }
    }
    // do the ghosts, on my level or finer
    for (auto biter = ghost_sibling_[qid->cid]->begin(); biter != ghost_sibling_[qid->cid]->end(); biter++) {
        GhostBlock* gblock = (*biter);
        // memory details
        MemLayout* block_src = ghost_subblock;
        real_p     data_src  = gblock->data_src();
        MemLayout* block_trg = gblock;
        real_p     data_trg  = cur_block->data(fid, ida_);
        // interpolate
        interp_->Interpolate(gblock->dlvl(), gblock->shift(), block_src, data_src, block_trg, data_trg);

        // we need to interpolate on the coarse version of myself as well
        if (do_coarse) {
            // set the coarse block to the correct position
            for (int id = 0; id < 3; id++) {
                coarse_start[id] = CoarseFromBlock(gblock->start(id));
                coarse_end[id]   = CoarseFromBlock(gblock->end(id));
            }
            coarse_subblock->Reset(M_GS, M_CLEN, coarse_start, coarse_end);
            // memory details
            MemLayout* block_src = ghost_subblock;
            real_p     data_src  = gblock->data_src();
            MemLayout* block_trg = coarse_subblock;
            real_p     data_trg  = tmp + m_zeroidx(0, coarse_subblock);
            // interpolate, the level is 1 coarser and the shift is unchanged
            interp_->Interpolate(gblock->dlvl() + 1, gblock->shift(), block_src, data_src, block_trg, data_trg);
        }
    }

    //-------------------------------------------------------------------------
    // copy the coarse blocks to the coarse representation
    for (auto biter = block_parent_[qid->cid]->begin(); biter != block_parent_[qid->cid]->end(); biter++) {
        m_verb("copying the coarse blocks");
        GhostBlock* gblock    = (*biter);
        GridBlock*  ngh_block = gblock->block_src();
        // setup the coarse sublock to the position
        for (int id = 0; id < 3; id++) {
            coarse_start[id] = CoarseFromBlock(gblock->start(id));
            coarse_end[id]   = CoarseFromBlock(gblock->end(id));
        }
        m_verb("originally: start = %d %d %d and end = %d %d %d", gblock->start(0), gblock->start(1), gblock->start(2), gblock->end(0), gblock->end(1), gblock->end(2));
        coarse_subblock->Reset(M_GS, M_CLEN, coarse_start, coarse_end);
        // memory details
        MemLayout* block_src = ngh_block;
        real_p     data_src  = ngh_block->data(fid, ida_);
        MemLayout* block_trg = coarse_subblock;
        real_p     data_trg  = tmp + m_zeroidx(0, coarse_subblock);
        // interpolate, the level is 1 coarser and the shift is unchanged
        m_assert(gblock->dlvl() + 1 == 0, "the gap in level has to be 0");
        interp_->Interpolate(gblock->dlvl() + 1, gblock->shift(), block_src, data_src, block_trg, data_trg);
    }
    // copy the ghost into the coarse representation
    for (auto biter = ghost_parent_[qid->cid]->begin(); biter != ghost_parent_[qid->cid]->end(); biter++) {
        m_verb("copying the coarse ghosts");
        GhostBlock* gblock = (*biter);
        // update the coarse subblock
        for (int id = 0; id < 3; id++) {
            coarse_start[id] = CoarseFromBlock(gblock->start(id));
            coarse_end[id]   = CoarseFromBlock(gblock->end(id));
        }
        coarse_subblock->Reset(M_GS, M_CLEN, coarse_start, coarse_end);
        // memory details
        MemLayout* block_src = ghost_subblock;
        real_p     data_src  = gblock->data_src();
        MemLayout* block_trg = coarse_subblock;
        real_p     data_trg  = tmp + m_zeroidx(0, coarse_subblock);
        // interpolate, the level is 1 coarser and the shift is unchanged
        m_assert(gblock->dlvl() + 1 == 0, "the gap in level has to be 0");
        interp_->Interpolate(gblock->dlvl() + 1, gblock->shift(), block_src, data_src, block_trg, data_trg);
    }
    //-------------------------------------------------------------------------
    // do a coarse version of myself and complete with some physics if needed
    if (do_coarse) {
        // set the coarse block to its whole domain
        for (int id = 0; id < 3; id++) {
            coarse_start[id] = 0;
            coarse_end[id]   = M_HN;
        }
        coarse_subblock->Reset(M_GS, M_CLEN, coarse_start, coarse_end);
        // get memory details
        lid_t      shift[3]  = {0, 0, 0};
        MemLayout* block_src = cur_block;
        real_p     data_src  = cur_block->data(fid, ida_);
        MemLayout* block_trg = coarse_subblock;
        real_p     data_trg  = tmp + m_zeroidx(0, coarse_subblock);
        // interpolate
        interp_->Interpolate(1, shift, block_src, data_src, block_trg, data_trg);

        // do here some physics

        // reset the coarse sublock to the full position
        for (int id = 0; id < 3; id++) {
            coarse_start[id] = -M_GS;
            coarse_end[id]   = M_HN + M_GS;
        }
        coarse_subblock->Reset(M_GS, M_CLEN, coarse_start, coarse_end);
    }
    //-------------------------------------------------------------------------
    // interpolate the ghost representation on myself
    for (auto biter = block_parent_[qid->cid]->begin(); biter != block_parent_[qid->cid]->end(); biter++) {
        GhostBlock* gblock    = (*biter);
        lid_t       shift[3]  = {0, 0, 0};
        MemLayout*  block_src = coarse_subblock;
        real_p      data_src  = tmp + m_zeroidx(0, coarse_subblock);
        MemLayout*  block_trg = gblock;
        real_p      data_trg  = cur_block->data(fid, ida_);
        // interpolate
        interp_->Interpolate(-1, shift, block_src, data_src, block_trg, data_trg);
    }
    // copy the ghost into the coarse representation
    for (auto biter = ghost_parent_[qid->cid]->begin(); biter != ghost_parent_[qid->cid]->end(); biter++) {
        GhostBlock* gblock    = (*biter);
        lid_t       shift[3]  = {0, 0, 0};
        MemLayout*  block_src = coarse_subblock;
        real_p      data_src  = tmp + m_zeroidx(0, coarse_subblock);
        MemLayout*  block_trg = gblock;
        real_p      data_trg  = cur_block->data(fid, ida_);
        // interpolate
        interp_->Interpolate(-1, shift, block_src, data_src, block_trg, data_trg);
    }

    //-------------------------------------------------------------------------
    // finally do some physics
    for (auto piter = phys_[qid->cid]->begin(); piter != phys_[qid->cid]->end(); piter++) {
        PhysBlock* gblock = (*piter);
        // get the direction and the corresponding bctype
        bctype_t bctype = fid->bctype(ida_, gblock->iface());
        if (bctype == M_BC_EVEN) {
            EvenBoundary_4 bc = EvenBoundary_4();
            bc(0.0, cur_block->hgrid(), gblock, cur_block->data(fid, ida_));
        } else if (bctype == M_BC_ODD) {
            OddBoundary_4 bc = OddBoundary_4();
            bc(0.0, cur_block->hgrid(), gblock, cur_block->data(fid, ida_));
        } else if (bctype == M_BC_EXTRAP) {
            ExtrapBoundary_4 bc = ExtrapBoundary_4();
            bc(0.0, cur_block->hgrid(), gblock, cur_block->data(fid, ida_));
        } else if (bctype == M_BC_ZERO) {
            ZeroBoundary bc = ZeroBoundary();
            bc(0.0, cur_block->hgrid(), gblock, cur_block->data(fid, ida_));
        } else {
            m_assert(false, "this type of BC is not implemented yet");
        }
    }
    delete (ghost_subblock);
    //-------------------------------------------------------------------------
}

void Ghost::LoopOnMirrorBlock_(const gop_t op, Field* field) {
    m_begin;
    m_assert(grid_->is_mesh_valid(), "mesh is not valid, unable to process");
    //-------------------------------------------------------------------------
    // get the grid info
    p8est_t*       forest = grid_->forest();
    p8est_ghost_t* ghost  = grid_->ghost();
    // const lid_t    nqlocal = ghost->mirrors.elem_count;  //number of ghost blocks

#pragma omp parallel for
    for (lid_t bid = 0; bid < n_mirror_to_send_; bid++) {
        // get the mirror quad, this is an empty quad (just a piggy3 struct)
        p8est_quadrant_t* mirror = p8est_quadrant_array_index(&ghost->mirrors, mirrors_to_local_[bid]);
        p8est_tree_t*     tree   = p8est_tree_array_index(forest->trees, mirror->p.piggy3.which_tree);
        // get the id, in this case the cummulative id = mirror id
        qid_t myid;
        myid.cid = bid;
        myid.qid = mirror->p.piggy3.local_num - tree->quadrants_offset;
        myid.tid = mirror->p.piggy3.which_tree;
        // use it to retreive the actual quadrant in the correct tree
        p8est_quadrant_t* quad  = p8est_quadrant_array_index(&tree->quadrants, myid.qid);
        GridBlock*        block = reinterpret_cast<GridBlock*>(quad->p.user_data);
        // send the task
        (this->*op)(&myid, block, field);
    }
    //-------------------------------------------------------------------------
    m_end;
}
