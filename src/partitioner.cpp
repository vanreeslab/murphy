#include "partitioner.hpp"

#include <cstring>
#include <map>

/**
 * @brief returns the rank that owns the value in the array of length commsize+1
 * 
 * @param array the array in which we search the value @ref val (length = commsize+1)
 * @param val the value we need to find in @ref array
 * @param range the search length
 * @param rankoffset where in the original array we started to look
 * @return int the position of @ref val
 */
static int bsearch_comm(const p4est_gloidx_t *array, const p4est_gloidx_t val, const int range, const int rankoffset) {
    const int middle = range / 2;
    if (array[middle] <= val && val < array[middle + 1]) {
        // we found it, return the corresponding rank
        return rankoffset + middle;
    } else if (val < array[middle]) {
        // we go for left
        // the new range is middle
        return bsearch_comm(array, val, middle, rankoffset);
    } else {
        // we go for right
        // the new range is range - middle
        return bsearch_comm(array + middle, val, range - middle, rankoffset + middle);
    }
}

using std::array;
using std::map;
using std::memcpy;

Partitioner::Partitioner(map<string, Field*>* fields, ForestGrid *grid) {
    m_begin;
    //-------------------------------------------------------------------------
    const int     commsize = grid->mpisize();
    p8est_t *     forest   = grid->forest();
    const int     rank     = forest->mpirank;

    // count how many ldas in total we have
    for (auto fid = fields->begin(); fid != fields->end(); fid++) {
        n_lda_ += fid->second->lda();
    }

    // store the lcoation in the old partition
    // note: we have to know the new partition to use it
    p4est_gloidx_t *oldpart = (p4est_gloidx_t *)m_calloc((commsize + 1) * sizeof(p4est_gloidx_t));
    memcpy(oldpart, forest->global_first_quadrant, (commsize + 1) * sizeof(p4est_gloidx_t));

    // init the array of current blocks
    const lid_t n_loc_block = forest->local_num_quadrants;
    old_blocks_             = (GridBlock **)m_calloc(n_loc_block * sizeof(GridBlock *));

    // store all the old block adresses before they get send
    for (p4est_topidx_t it = forest->first_local_tree; it <= forest->last_local_tree; it++) {
        p8est_tree_t *mytree = p8est_tree_array_index(forest->trees, it);
        for (lid_t qid = 0; qid < mytree->quadrants.elem_count; qid++) {
            qdrt_t *quad   = p8est_quadrant_array_index(&mytree->quadrants, qid);
            lid_t   offset = mytree->quadrants_offset;
            // store the block adress
            old_blocks_[offset + qid] = *(reinterpret_cast<GridBlock**>(quad->p.user_data));
        }
    }

    // compute the new partition, asking the children to be on the same block (in case of coarsening)
    p4est_gloidx_t nqshipped = p8est_partition_ext(forest, 0, nullptr);

    m_log("total %ld blocks are moving: going from %d to now %d in here",nqshipped,n_loc_block,forest->local_num_quadrants);

    if (nqshipped > 0) {
        // get the NEW number of quads
        const lid_t nqlocal = forest->local_num_quadrants;
        // get the new starting point and end point of the curent partition
        const p4est_gloidx_t opart_begin = oldpart[rank];
        const p4est_gloidx_t opart_end   = oldpart[rank + 1];
        const p4est_gloidx_t cpart_begin = forest->global_first_quadrant[rank];
        const p4est_gloidx_t cpart_end   = forest->global_first_quadrant[rank + 1];
        // count how many blocks are common, so that they don't travel
        const lid_t q_nself = m_max(0, m_min(opart_end, cpart_end) - m_max(opart_begin, cpart_begin));
        const lid_t opart_n = m_max(0, opart_end - opart_begin - q_nself);
        const lid_t cpart_n = m_max(0, cpart_end - cpart_begin - q_nself);

        m_log("except the self, I lose %d blocks and gain %d blocks",opart_n,cpart_n);

        // init the send
        if (opart_n > 0) {
            // the send buffer is used to copy the current blocks as they are not continuous to memory
            send_buf_ = (real_t *)m_calloc(opart_n * m_blockmemsize(n_lda_) * sizeof(real_t));

            //receivers = from opart_begin to opart_end-1 in the new partition
            const int first_recver = bsearch_comm(forest->global_first_quadrant, opart_begin, commsize, 0);
            const int last_recver  = bsearch_comm(forest->global_first_quadrant, opart_end - 1, commsize, 0);
            const int n_recver     = last_recver - first_recver + 1;  // note:add one because 4-2=2 but we have 3 ranks...
            // if I am among the list, remove myself from the request point of view
            n_send_request_ = (first_recver <= rank && rank <= last_recver) ? n_recver - 1 : n_recver;
            // allocate the requests
            send_request_ = (MPI_Request *)m_calloc(n_send_request_ * sizeof(MPI_Request));
            // allocate the cummulative list, to know which block I have to copy to the buffer
            // note: we cannot use a cummulative since we may not have continuous set (if some blocks are not moving)
            q_send_start_block_ = (lid_t *)m_calloc(n_send_request_ * sizeof(lid_t));
            q_send_cum_request_   = (lid_t *)m_calloc((n_send_request_ + 1) * sizeof(lid_t));

            // for each receiver, allocate the send request
            int   rcount = 0;
            lid_t qcount = 0;
            lid_t tqcount = 0;
            for (int ir = 0; ir < n_recver; ir++) {
                // get who is the receiver and skip if it's me
                const int c_recver = first_recver + ir;
                // get howmany the receiver will receive (not his/her entire block numbers)
                const p4est_gloidx_t q_leftlimit  = m_max(forest->global_first_quadrant[c_recver], opart_begin);
                const p4est_gloidx_t q_rightlimit = m_min(forest->global_first_quadrant[c_recver + 1], opart_end);
                const lid_t          n_q2send     = q_rightlimit - q_leftlimit;
                // if I want to allocate a send for myself, skip
                if (c_recver == rank) {
                    tqcount += n_q2send;
                    continue;
                }
                // remember the begin and end point
                q_send_start_block_[rcount] = tqcount;
                q_send_cum_request_[rcount] = qcount + n_q2send;

                // create the send request
                real_t *buf = send_buf_ + qcount * m_blockmemsize(n_lda_);
                MPI_Send_init(buf, m_blockmemsize(n_lda_) * n_q2send, M_MPI_REAL, c_recver, c_recver, forest->mpicomm, &(send_request_[rcount]));

                // update the counters
                qcount  += n_q2send;
                tqcount += n_q2send;
                rcount  += 1;
            }
            // store the last index
            q_send_cum_request_[n_send_request_] = qcount;
            m_assert(qcount == opart_n, "counters are not matching");

            m_log("done with %d send requests",n_send_request_);
        }

        // init the reception
        if (cpart_n > 0) {
            // store the new quadrant adress, to create a new block if needed
            new_blocks_  = (GridBlock **)m_calloc(forest->local_num_quadrants * sizeof(GridBlock *));
            p4est_gloidx_t rank_offset = forest->global_first_quadrant[rank];

            // got through each block and add it if needed
            lid_t bcount = 0;
            for (p4est_topidx_t it = forest->first_local_tree; it <= forest->last_local_tree; it++) {
                p8est_tree_t *mytree = p8est_tree_array_index(forest->trees, it);
                for (lid_t qid = 0; qid < mytree->quadrants.elem_count; qid++) {
                    qdrt_t *quad   = p8est_quadrant_array_index(&mytree->quadrants, qid);
                    lid_t   offset = mytree->quadrants_offset;
                    // get the global ID and allocate a new block only if the quad is new
                    p4est_gloidx_t global_id = rank_offset + offset + qid;
                    if (!(opart_begin <= global_id && global_id < opart_end)) {
                        // create a new block
                        real_t xyz[3];
                        p8est_qcoord_to_vertex(forest->connectivity, it, quad->x, quad->y, quad->z, xyz);
                        real_t     len   = m_quad_len(quad->level);
                        GridBlock *block = new GridBlock(len, xyz, quad->level);
                        // add the fields to the block
                        block->AddFields(fields);
                        // store its access, replace the adress that was there but which is wrong now
                        *(reinterpret_cast<GridBlock **>(quad->p.user_data)) = block;
                        // store in the array
                        new_blocks_[offset + qid] = block;
                        // increment the counter
                        bcount ++;
                    } else {
                        new_blocks_[offset + qid] = nullptr;// *(reinterpret_cast<GridBlock **>(quad->p.user_data));
                    }
                }
            }
            m_assert(bcount == cpart_n,"the counters should match");

            // the receive buffer is used as a new data location for the blocks
            recv_buf_ = (real_t *)m_calloc(cpart_n * m_blockmemsize(n_lda_) * sizeof(real_t));

            //senders = from cpart_begin to cpart_end-1 in the new partition
            m_log("looking for the senders of %ld -> %ld",cpart_begin,cpart_end);
            const int first_sender = bsearch_comm(oldpart, cpart_begin, commsize, 0);
            const int last_sender  = bsearch_comm(oldpart, cpart_end - 1, commsize, 0);
            m_log("found sender %d to %d",first_sender,last_sender);
            const int n_sender     = last_sender - first_sender + 1;  // note:add one because 4-2=2 but we have 3 ranks...
            // if I am among the list, remove myself from the request point of view
            n_recv_request_ = (first_sender <= rank && rank <= last_sender) ? n_sender - 1 : n_sender;
            m_log("I will do %d recv reqests our of %d senders", n_recv_request_, n_sender);
            // allocate the requests
            recv_request_ = (MPI_Request *)m_calloc(n_recv_request_ * sizeof(MPI_Request));
            // allocate the cummulative list, to know which block I have to copy to the buffer
            // note: we cannot use a cummulative since we may not have continuous set (if some blocks are not moving)
            q_recv_start_block_ = (lid_t *)m_calloc(n_recv_request_ * sizeof(lid_t));
            q_recv_cum_request_ = (lid_t *)m_calloc((n_recv_request_ + 1) * sizeof(lid_t));

            // for each receiver, allocate the send request
            int   scount  = 0;
            lid_t qcount  = 0;  //! counter on the blocks actually send
            lid_t tqcount = 0;  //! counter on all the blocks on this rank
            for (int ir = 0; ir < n_sender; ir++) {
                // get who is the sender
                const int c_sender = first_sender + ir;

                // get howmany the sender will send (not his/her entire block numbers)
                const p4est_gloidx_t q_leftlimit  = m_max(oldpart[c_sender], cpart_begin);
                const p4est_gloidx_t q_rightlimit = m_min(oldpart[c_sender + 1], cpart_end);
                const lid_t          n_q2recv     = q_rightlimit - q_leftlimit;
                // if I am the sender, advance the tqcount counter and skip
                if (c_sender == rank) {
                    tqcount += n_q2recv;
                    continue;
                }
                // store the memory accesses
                q_recv_start_block_[scount] = tqcount;
                q_recv_cum_request_[scount] = qcount + n_q2recv;
                // create the send request
                real_p buf = recv_buf_ + qcount * m_blockmemsize(n_lda_);
                MPI_Recv_init(buf, m_blockmemsize(n_lda_) * n_q2recv, M_MPI_REAL, c_sender, rank, forest->mpicomm, &(recv_request_[scount]));
                // update the counters
                qcount += n_q2recv;
                tqcount += n_q2recv;
                scount += 1;
            }
            q_recv_cum_request_[n_recv_request_] = qcount;
            m_assert(qcount == cpart_n, "counters are not matching");
        }
    }
    m_free(oldpart);

    m_log("grid partitioned: moving %ld blocks (%f percent of the grid): %d sends and %d recvs",nqshipped,((real_t)nqshipped)/((real_t)forest->global_num_quadrants),n_send_request_,n_recv_request_);

    //-------------------------------------------------------------------------
    m_end;
}

Partitioner::~Partitioner() {
    DeallocOldies_();
    if (old_blocks_ != nullptr) {
        m_free(old_blocks_);
    }
    if (new_blocks_ != nullptr) {
        m_free(new_blocks_);
    }

    for (int i = 0; i < n_send_request_; i++) {
        MPI_Request_free(send_request_ + i);
    }
    for (int i = 0; i < n_recv_request_; i++) {
        MPI_Request_free(recv_request_ + i);
    }

    m_free(send_request_);
    m_free(recv_request_);
    m_free(send_buf_);
    m_free(recv_buf_);
    m_free(q_send_start_block_);
    m_free(q_send_cum_request_);
    m_free(q_recv_start_block_);
    m_free(q_recv_cum_request_);
}

void Partitioner::Start(map<string, Field *> *fields) {
    m_begin;
    //-------------------------------------------------------------------------
    // start the reception of the memory
    if (n_recv_request_ > 0) {
        MPI_Startall(n_recv_request_, recv_request_);
    }
    // copy all the data of the old blocks into the send_buffer and send when ready
    for (int is = 0; is < n_send_request_; is++) {
        lid_t n_q2send = q_send_cum_request_[is + 1] - q_send_cum_request_[is];
        for (lid_t iq = 0; iq < n_q2send; iq++) {
            GridBlock *block = old_blocks_[q_send_start_block_[is] + iq];
            real_p     buf   = send_buf_ + (q_send_cum_request_[is] + iq) * m_blockmemsize(n_lda_);

            lid_t idacount = 0;
            for (auto iter = fields->begin(); iter != fields->end(); fields++) {
                const Field *fid  = iter->second;
                string       name = fid->name();
                // the memory returnded by block->data is shifted, so we unshift it :-)
                real_p data = block->data(fid) - m_zeroidx(0, block);
                real_p lbuf = buf + m_blockmemsize(idacount);
                // copy the whole memory at once on the dim
                memcpy(lbuf, data, m_blockmemsize(fid->lda()) * sizeof(real_t));
                // update the counters
                idacount += fid->lda();
            }
        }
        // start the send
        MPI_Start(send_request_ + is);
    }
    //-------------------------------------------------------------------------
    m_end;
}

void Partitioner::End(map<string, Field *> *fields) {
    for (int ir = 0; ir < n_recv_request_; ir++) {
        int idx;
        MPI_Waitany(n_recv_request_, recv_request_, &idx, MPI_STATUS_IGNORE);
        lid_t n_q2recv = q_recv_cum_request_[idx + 1] - q_recv_cum_request_[idx];

        for (lid_t iq = 0; iq < n_q2recv; iq++) {
            GridBlock *block = new_blocks_[q_recv_start_block_[idx] + iq];
            real_p     buf   = recv_buf_ + (q_recv_cum_request_[idx] + iq) * m_blockmemsize(n_lda_);
            m_assert(block != nullptr, "this block shouldn't be accessed here");

            lid_t idacount = 0;
            for (auto iter = fields->begin(); iter != fields->end(); fields++) {
                const Field *fid  = iter->second;
                string       name = fid->name();
                // the memory returnded by block->data is shifted, so we unshift it :-)
                real_p data = block->data(fid) - m_zeroidx(0, block);
                real_p lbuf = buf + m_blockmemsize(idacount);
                // copy the whole memory at once on the dim
                memcpy(data, lbuf, m_blockmemsize(fid->lda()) * sizeof(real_t));
                // update the counters
                idacount += fid->lda();
            }
        }
    }
    // receive the new memory in the new blocks
    if (n_send_request_ > 0) {
        MPI_Waitall(n_send_request_, send_request_, MPI_STATUS_IGNORE);
    }
}

void Partitioner::DeallocOldies_() {
    m_begin;
    //-------------------------------------------------------------------------
    for (int is = 0; is < n_send_request_; is++) {
        // delete all the blocks that have been transfered, not the remaining ones..
        lid_t n_q2send = q_send_cum_request_[is + 1] - q_send_cum_request_[is];
        for (lid_t iq = q_send_start_block_[is]; iq < (q_send_start_block_[is] + n_q2send); iq++) {
            if (old_blocks_[iq] != nullptr) {
                delete (old_blocks_[iq]);
            }
        }
    }
    //-------------------------------------------------------------------------
    m_end;
}