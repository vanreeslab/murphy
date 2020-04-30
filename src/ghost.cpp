#include "ghost.hpp"

#include <limits>

#include "boundary.hpp"
#include "murphy.hpp"
#include "omp.h"
#include "wavelet.hpp"

/**
 * @brief returns the number of ghost points for the coarse block
 * 
 * The ghost are padded so that the array is aligned
 */
static lid_t cgs(Interpolator* interp) {
    lid_t nghost    = (interp->NGhostCoarse() + 1);
    lid_t n_in_line = (M_ALIGNMENT / sizeof(real_t));
    nghost          = nghost + n_in_line - (nghost % n_in_line);
    return nghost;
}
/**
 * @brief returns the stride of the coarse block
 * it needs NGhostCoarse()+1 ghost points to perform the refinement in the fine's ghost points
 */
static size_t cstride(Interpolator* interp) {
    return 2 * cgs(interp) + M_HN;
}

/**
 * @brief Construct a new Ghost given a ForestGrid and an interpolator, allocate the ghost lists and initiates the communications
 * 
 * Once created, the ghost is fixed for a given grid. if the grid changes, a new Ghost objects has to be created.
 * 
 * @note: The Interpolator has to be given beforehands because the output of Interpolator::NGhostCoarse() drives the coarse block allocation
 * 
 * @param grid
 * @param interp the interpolator used to interpolate the ghost values
 */
Ghost::Ghost(ForestGrid* grid, Interpolator* interp) {
    m_begin;
    m_assert(grid->is_mesh_valid(), "the mesh needs to be valid before entering here");
    //-------------------------------------------------------------------------
    // get the important pointers
    grid_              = grid;
    interp_            = interp;
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
    coarse_tmp_  = reinterpret_cast<real_t**>(m_calloc(sizeof(real_t*) * nthreads));
    for (int it = 0; it < nthreads; it++) {
        coarse_tmp_[it] = reinterpret_cast<real_t*>(m_calloc(sizeof(real_t) * cstride(interp_) * cstride(interp_) * cstride(interp_)));
    }
    //-------------------------------------------------------------------------
    // m_log("ghost initialized with interpolator: %s -> %f percent of ghost comm not needed",interp_->Identity().c_str(),100.0*pow((M_N-2*cgs(interp_))/((real_t)M_N),3));
    m_log("ghost initialized with %s", interp_->Identity().c_str());
    m_end;
}

/**
 * @brief Destroy the Ghost and free the allocated memory
 * 
 */
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
    n_send_request_ = 0;
    for (int ir = 0; ir < n_recv_request_; ir++) {
        MPI_Request_free(ghost_recv_ + ir);
    }
    m_free(ghost_recv_);
    ghost_recv_ = 0;

    // free the temp memory
    int nthreads = omp_get_max_threads();
    for (int it = 0; it < nthreads; it++) {
        m_free(coarse_tmp_[it]);
    }
    m_free(coarse_tmp_);

    // free the memory
    m_free(mirrors_);
    m_free(ghosts_);
    m_free(mirrors_to_local_);
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief given the dimension and the field, copy the mirror block's value to the send buffers
 * 
 * @param field 
 * @param ida 
 */
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

/**
 * @brief compute the ghost values from the received buffers and the local blocks
 * 
 * @param field 
 * @param ida 
 * @param interp 
 */
void Ghost::PullFromGhost(Field* field, sid_t ida) {
    m_begin;
    m_assert(grid_->is_mesh_valid(), "the mesh needs to be valid before entering here");
    //-------------------------------------------------------------------------
    // store the current interpolator and dimension
    ida_ = ida;
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
void Ghost::ApplyOpS(const qid_t* qid, GridBlock* block) {
    //-------------------------------------------------------------------------
    InitList_(qid, block);
    //-------------------------------------------------------------------------
}

/**
 * @brief Pull (= get the value for the neighbors) the ghost points for a block, given a field and the direction ida_
 * 
 * @param qid 
 * @param block 
 * @param fid 
 */
void Ghost::ApplyOpF(const qid_t* qid, GridBlock* cur_block, Field* fid) {
    //-------------------------------------------------------------------------
    PullFromGhost_(qid, cur_block, fid);
    //-------------------------------------------------------------------------
}

/**
 * @brief initalize the MPI_Request communication once, ready to be used
 * 
 * @note One p4est mirror structure can be send to multiple ranks
 * 
 */
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
    mirror_send_ = reinterpret_cast<MPI_Request*>(m_calloc(n_send_request_ * sizeof(MPI_Request)));
    ghost_recv_  = reinterpret_cast<MPI_Request*>(m_calloc(n_recv_request_ * sizeof(MPI_Request)));

    // allocate the mirror and ghost arrayss
    mirrors_to_local_ = reinterpret_cast<lid_t*>(m_calloc(sizeof(lid_t) * n_mirror_to_send_));
    mirrors_          = reinterpret_cast<real_t*>(m_calloc(sizeof(real_t) * M_NGHOST * n_mirror_to_send_));
    ghosts_           = reinterpret_cast<real_t*>(m_calloc(sizeof(real_t) * M_NGHOST * ghost->ghosts.elem_count));

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

/**
 * @brief setup the ghosting lists, for the considered block
 * 
 * @warning this function is called within a omp parellel region.
 * The p4est arrays are NOT thread-safe, while the lists have been assumed to be NOT threadsafe as well.
 * 
 * @note the p4est array allocation is thread-safe, while the tracking of the memory allocated is NOT.
 * 
 * @param qid the quadrant ID
 * @param block the associated block
 */
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
            // m_log("I have a grid with %d blocks and %d ghosts. Looking for neighbor %d of block %d = %d, %d",mesh->local_num_quadrants,ghost->ghosts.elem_count,ibidule,qid->cid,qid->qid,qid->tid);
            p8est_mesh_get_neighbors(forest, ghost, mesh, qid->cid, ibidule, ngh_quad, ngh_enc, ngh_qid);
        }
        // decode the status and count the ghosts
        const lid_t nghosts = ngh_enc->elem_count;
        m_assert(ngh_enc->elem_count < numeric_limits<lid_t>::max(), "the number of ghost is too big");
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
        for (lid_t nid = 0; nid < nghosts; nid++) {
            const int  status  = *(ngh_enc->array + nid * sizeof(int));
            const bool isghost = (status < 0);
            qdrt_t*    nghq    = *(reinterpret_cast<qdrt_t**>(sc_array_index_int(ngh_quad, nid)));

            // get the sign, i.e. the normal to the face, the edge of the corner we consider
            real_t sign[3];
            GhostGetSign(ibidule, sign);

            // get the position of the neighbor, as seen by me!!! may be different than the actual one if there is a periodic bc
            real_t ngh_pos[3];
            if (!isghost) {
                // cannot use the p8est function because the which_tree is not assigned, so we retrieve the position through the block
                GridBlock* ngh_block = *(reinterpret_cast<GridBlock**>(nghq->p.user_data));
                ngh_pos[0]           = ngh_block->xyz(0);
                ngh_pos[1]           = ngh_block->xyz(1);
                ngh_pos[2]           = ngh_block->xyz(2);
            } else {
                p8est_qcoord_to_vertex(grid_->connect(), nghq->p.piggy3.which_tree, nghq->x, nghq->y, nghq->z, ngh_pos);
            }
            for (sid_t id = 0; id < 3; id++) {
                ngh_pos[id] = block->xyz(id) + fmod(ngh_pos[id] - block->xyz(id) + grid_->domain_periodic(id) * sign[id] * grid_->domain_length(id), grid_->domain_length(id));
            }

            // create the new block
            GhostBlock* gb = new GhostBlock(block, nghq->level, ngh_pos);

            // associate the correct memory and push back
            if (!isghost) {
                GridBlock* ngh_block = *(reinterpret_cast<GridBlock**>(nghq->p.user_data));
                gb->block_src(ngh_block);

                if (gb->dlvl() >= 0) {
#pragma omp critical
                    bsibling->push_back(gb);
                } else {
#pragma omp critical
                    bparent->push_back(gb);
                }
            } else {
                lid_t ighost = *(reinterpret_cast<int*>(sc_array_index_int(ngh_qid, nid)));
                m_assert((ighost >= 0) && (ighost < ghost->ghosts.elem_count), "treeid = %d, qid = %d, ibidule = %d, the ID of the ghost is INVALID: %d vs %ld, status = %d (nid = %d, nghost = %d, array length=%ld)", qid->tid, qid->qid, ibidule, ighost, ghost->ghosts.elem_count, status, nid, nghosts, ngh_qid->elem_count);
                real_p data = ghosts_ + ighost * M_NGHOST;
                gb->data_src(data);

                if (gb->dlvl() >= 0) {
#pragma omp critical
                    gsibling->push_back(gb);
                } else {
#pragma omp critical
                    gparent->push_back(gb);
                }
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

/**
 * @brief copy the block information into the mirror buffer for a given field and dimension
 * 
 * @param qid the quadrant ID
 * @param block th block that is an actual mirror
 * @param fid the field ID
 */
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
 * @brief given an up-to-date @ref ghost_ arrays,
 * compute the ghosts points for the current block (field @ref fid, dimension @ref ida_)
 * 
 * @warning this function runs inside an OpenMP parallel region
 * 
 * @note a lot of dependencies appear if the block has a coarser neighbor. To solve it, we build
 * step by step a coarse representation of the current block, which is then used to refine to the actual block.
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
    SubBlock* coarse_subblock = new SubBlock(cgs(interp_), cstride(interp_) , coarse_start, coarse_end);

    // determine if we have to use the coarse representationun
    const bool do_coarse = (block_parent_[qid->cid]->size() + ghost_parent_[qid->cid]->size()) > 0;
    // if so, reset the coarse info
    if (do_coarse) {
        memset(tmp, 0, cstride(interp_) * cstride(interp_) * cstride(interp_) * sizeof(real_t));
    }

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
            coarse_subblock->Reset(cgs(interp_), cstride(interp_), coarse_start, coarse_end);
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
            coarse_subblock->Reset(cgs(interp_), cstride(interp_), coarse_start, coarse_end);
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
        GhostBlock* gblock    = (*biter);
        GridBlock*  ngh_block = gblock->block_src();
        // setup the coarse sublock to the position
        for (int id = 0; id < 3; id++) {
            coarse_start[id] = CoarseFromBlock(gblock->start(id));
            coarse_end[id]   = CoarseFromBlock(gblock->end(id));
        }
        coarse_subblock->Reset(cgs(interp_), cstride(interp_), coarse_start, coarse_end);
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
        GhostBlock* gblock = (*biter);
        // update the coarse subblock
        for (int id = 0; id < 3; id++) {
            coarse_start[id] = CoarseFromBlock(gblock->start(id));
            coarse_end[id]   = CoarseFromBlock(gblock->end(id));
        }
        coarse_subblock->Reset(cgs(interp_), cstride(interp_), coarse_start, coarse_end);
        // memory details
        MemLayout* block_src = ghost_subblock;
        real_p     data_src  = gblock->data_src();
        MemLayout* block_trg = coarse_subblock;
        real_p     data_trg  = tmp + m_zeroidx(0, coarse_subblock);
        // interpolate, the level is 1 coarser and the shift is unchanged
        m_assert(gblock->dlvl() + 1 == 0, "the gap in level has to be 0");
        interp_->Interpolate(gblock->dlvl() + 1, gblock->shift(), block_src, data_src, block_trg, data_trg);
    }

    // m_log("dbg: tree %d, quad %d: doing coarse ? %d", qid->tid, qid->qid,do_coarse);
    //-------------------------------------------------------------------------
    // do a coarse version of myself and complete with some physics if needed
    if (do_coarse) {
        // set the coarse block to its whole domain
        for (int id = 0; id < 3; id++) {
            coarse_start[id] = 0;
            coarse_end[id]   = M_HN;
        }
        coarse_subblock->Reset(cgs(interp_), cstride(interp_), coarse_start, coarse_end);
        // get memory details
        lid_t      shift[3]  = {0, 0, 0};
        MemLayout* block_src = cur_block;
        real_p     data_src  = cur_block->data(fid, ida_);
        MemLayout* block_trg = coarse_subblock;
        real_p     data_trg  = tmp + m_zeroidx(0, coarse_subblock);
        // interpolate
        interp_->Interpolate(1, shift, block_src, data_src, block_trg, data_trg);

        // do here some physics, to completely fill the coarse block before the interpolation
        for (auto piter = phys_[qid->cid]->begin(); piter != phys_[qid->cid]->end(); piter++) {
            PhysBlock* gblock = (*piter);
            // get the direction and the corresponding bctype
            const sid_t    dir    = gblock->dir();
            const bctype_t bctype = fid->bctype(ida_, gblock->iface());
            // in the face direction, the start and the end are already correct, only the fstart changes
            lid_t fstart[3];
            coarse_start[dir] = gblock->start(dir);
            coarse_end[dir]   = gblock->end(dir);
            fstart[dir]       = CoarseFromBlock(face_start[gblock->iface()][dir]);
            // in the other direction, we need to rescale the dimensions
            for (int id = 1; id < 3; id++) {
                coarse_start[(dir + id) % 3] = CoarseFromBlock(gblock->start((dir + id) % 3));
                coarse_end[(dir + id) % 3]   = CoarseFromBlock(gblock->end((dir + id) % 3));
                fstart[(dir + id) % 3]       = CoarseFromBlock(face_start[gblock->iface()][(dir + id) % 3]);
            }
            // reset the coarse block and get the correct memory location
            coarse_subblock->Reset(cgs(interp_), cstride(interp_), coarse_start, coarse_end);
            real_p data_trg = tmp + m_zeroidx(0, coarse_subblock);
            // get the correct face_start
            if (bctype == M_BC_EVEN) {
                EvenBoundary_4 bc = EvenBoundary_4();
                bc(gblock->iface(), fstart, cur_block->hgrid(), 0.0, coarse_subblock, data_trg);
            } else if (bctype == M_BC_ODD) {
                OddBoundary_4 bc = OddBoundary_4();
                bc(gblock->iface(), fstart, cur_block->hgrid(), 0.0, coarse_subblock, data_trg);
            } else if (bctype == M_BC_EXTRAP) {
                ExtrapBoundary_4 bc = ExtrapBoundary_4();
                bc(gblock->iface(), fstart, cur_block->hgrid(), 0.0, coarse_subblock, data_trg);
            } else if (bctype == M_BC_ZERO) {
                ZeroBoundary bc = ZeroBoundary();
                bc(gblock->iface(), fstart, cur_block->hgrid(), 0.0, coarse_subblock, data_trg);
            } else {
                m_assert(false, "this type of BC is not implemented yet");
            }
        }

        // reset the coarse sublock to the full position
        for (int id = 0; id < 3; id++) {
            coarse_start[id] = -cgs(interp_);
            coarse_end[id]   = cstride(interp_)-cgs(interp_);
        }
        coarse_subblock->Reset(cgs(interp_), cstride(interp_), coarse_start, coarse_end);
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
            bc(gblock->iface(), face_start[gblock->iface()], cur_block->hgrid(), 0.0, gblock, cur_block->data(fid, ida_));
        } else if (bctype == M_BC_ODD) {
            OddBoundary_4 bc = OddBoundary_4();
            bc(gblock->iface(), face_start[gblock->iface()], cur_block->hgrid(), 0.0, gblock, cur_block->data(fid, ida_));
        } else if (bctype == M_BC_EXTRAP) {
            ExtrapBoundary_4 bc = ExtrapBoundary_4();
            bc(gblock->iface(), face_start[gblock->iface()], cur_block->hgrid(), 0.0, gblock, cur_block->data(fid, ida_));
        } else if (bctype == M_BC_ZERO) {
            ZeroBoundary bc = ZeroBoundary();
            bc(gblock->iface(), face_start[gblock->iface()], cur_block->hgrid(), 0.0, gblock, cur_block->data(fid, ida_));
        } else {
            m_assert(false, "this type of BC is not implemented yet");
        }
    }
    delete (ghost_subblock);
    delete (coarse_subblock);

    // m_log("dbg: tree %d, quad %d: finiiiish", qid->tid, qid->qid);
    //-------------------------------------------------------------------------
}

/**
 * @brief Loop on the blocks that are mirrors and call a gop_t operation on them
 * 
 * @param op 
 * @param field 
 */
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
        GridBlock*        block = *(reinterpret_cast<GridBlock**>(quad->p.user_data));
        // send the task
        (this->*op)(&myid, block, field);
    }
    //-------------------------------------------------------------------------
    m_end;
}
