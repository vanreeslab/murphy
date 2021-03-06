#include "partitioner.hpp"

#include <cstring>
#include <map>
#include <unordered_map>

using std::string;
using std::unordered_map;

/**
 * @brief returns the rank that owns the value in the array of length commsize+1
 * 
 * @param array the array in which we search the value @ref val (length = commsize+1)
 * @param val the value we need to find in @ref array
 * @param range the search length
 * @param rankoffset where in the original array we started to look
 * @return int the position of @ref val
 */
static rank_t bsearch_comm(const p4est_gloidx_t *array, const p4est_gloidx_t val, const int range, const int rankoffset) {
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

#ifndef NDEBUG
static size_t PartitionerBlockSize(const CartBlock *block, const lda_t lda) {
    const size_t n_elem = block->BlockLayout().n_elem;
    const size_t offset = block->PartitionDataOffset();
    return (n_elem * lda + offset);
}
#endif

/**
 * @brief Construct a new Partitioner by changing the given ForestGrid to match the new partition and register the communication pattern
 * 
 * @warning after the execution of it: the grid follows the new parition, new, yet empty, GridBlock have been allocated on the new partition and
 * the information is still on the old GridBlock on the old partition. Hence a significant memory overhead can be observed.
 * 
 * @param fields the list of fields to partition, we reallocate the only that list of the fields
 * @param grid the grid on which the partitionning is done
 * @param destructive dictate if the partitioning is used to permanently change the grid or to go from one grid to another.
 */
Partitioner::Partitioner(map<string, Field *> *fields, Grid *grid, bool destructive) {
    m_begin;
    //-------------------------------------------------------------------------
    const rank_t commsize = grid->mpisize();
    p8est_t *    forest   = grid->p4est_forest();
    const rank_t rank     = forest->mpirank;
    m_assert(commsize >= 0, "the commsize = %d must be positive", commsize);

    // store the destructive state
    destructive_ = destructive;

    // count how many ldas we have to forecast for the exchange, depend on the destructive mode
    n_lda_ = 0;
    if (destructive_) {
        // if destructive, we move all the fields
        m_assert(fields->size() == grid->NField(), "the number of fields must match so we don't loose information during the partitioning");
        for (const auto fid : *fields) {
            // add the total count
            n_lda_ += fid.second->lda() * (!fid.second->is_expr());
            m_assert(grid->IsAField(fid.second), "the field MUST be present in the grid");
        }
    } else {
        // create the struct to send the max out of the fields
        for (const auto fid : *fields) {
            n_lda_ = m_max(n_lda_, fid.second->lda() * (!fid.second->is_expr()));
            // m_assert(grid->IsAField(fid->second),"the field MUST be present in the grid");
        }
    }
    // we might have to send no field, but never a negative one!!
    m_assert(n_lda_ >= 0, "I cannot send a negative dimension");

    //..........................................................................
    // store the location of the quads in the old partition
    // note: we have to know the new partition to use it

    p4est_gloidx_t *oldpart = reinterpret_cast<p4est_gloidx_t *>(m_calloc((commsize + 1) * sizeof(p4est_gloidx_t)));
    memcpy(oldpart, forest->global_first_quadrant, (commsize + 1) * sizeof(p4est_gloidx_t));

    // get the blocksize from the grid. We cannot get it from the quadrants or the blocks
    // because we might have no quads and still get some coming.
    partition_blocksize_ = BlockPartitionSize(grid->block_type(), n_lda_);

    // init the array of current blocks and store their adress before we send it
    const iblock_t n_loc_block = forest->local_num_quadrants;
    old_blocks_                = reinterpret_cast<CartBlock **>(m_calloc(n_loc_block * sizeof(CartBlock *)));
    // store all the old block adresses before they get send
    for (p4est_topidx_t it = forest->first_local_tree; it <= forest->last_local_tree; ++it) {
        p8est_tree_t *mytree = p8est_tree_array_index(forest->trees, it);
        for (iblock_t qid = 0; qid < mytree->quadrants.elem_count; ++qid) {
            const qdrt_t * quad   = p8est_quadrant_array_index(&mytree->quadrants, qid);
            const iblock_t offset = mytree->quadrants_offset;
            // store the block address
            CartBlock *block          = p4est_GetBlock<CartBlock>(quad);
            old_blocks_[offset + qid] = block;
            m_assert(partition_blocksize_ == PartitionerBlockSize(block, n_lda_), "The block size must be constant over the grid: %ld vs %ld", partition_blocksize_, PartitionerBlockSize(block, n_lda_));
        }
    }

    //..........................................................................
    m_verb("current status: %d quads locally", forest->local_num_quadrants);
    // compute the new partition, asking the children to be on the same block (in case of coarsening)
    p4est_gloidx_t nqshipped = p8est_partition_ext(forest, true, NULL);
    m_assert(nqshipped >= 0, "the number of quads to send must be >= 0");
    m_verb("we decided to move %ld blocks", nqshipped);
    m_verb("new status: %d quads locally", forest->local_num_quadrants);

    //..........................................................................
    if (nqshipped > 0) {
        // get the NEW number of quads
        const iblock_t nqlocal = forest->local_num_quadrants;
        // get the new starting point and end point of the curent partition
        const p4est_gloidx_t opart_begin = oldpart[rank];
        const p4est_gloidx_t opart_end   = oldpart[rank + 1];
        const p4est_gloidx_t cpart_begin = forest->global_first_quadrant[rank];
        const p4est_gloidx_t cpart_end   = forest->global_first_quadrant[rank + 1];
        // count how many blocks are common, so that they don't travel
        m_verb("I had before %ld to %ld and now %ld to %ld", opart_begin, opart_end, cpart_begin, cpart_end);
        const iblock_t q_nself = m_max(0, m_min(opart_end, cpart_end) - m_max(opart_begin, cpart_begin));
        const iblock_t opart_n = m_max(0, opart_end - opart_begin - q_nself);
        const iblock_t cpart_n = m_max(0, cpart_end - cpart_begin - q_nself);

        m_verb("except the self (= %d blocks), I lose %d blocks and gain %d blocks", q_nself, opart_n, cpart_n);

        //......................................................................
        // init the send
        if (opart_n > 0) {
            // the send buffer is used to copy the current blocks as they are not continuous to memory
            send_buf_.Allocate(opart_n * partition_blocksize_);
            // if (destructive_) {
            //     // the status are only send if the partitioner is destructive
            //     send_status_buf_ = reinterpret_cast<short *>(m_calloc(opart_n * sizeof(short)));
            // }
            m_verb("sending buffer initialize of size %ld bytes", opart_n * partition_blocksize_ * sizeof(real_t));

            // receivers = from opart_begin to opart_end-1 in the new partition
            const rank_t first_recver = bsearch_comm(forest->global_first_quadrant, opart_begin, commsize, 0);
            const rank_t last_recver  = bsearch_comm(forest->global_first_quadrant, opart_end - 1, commsize, 0);
            const rank_t n_recver     = last_recver - first_recver + 1;  // note:add one because 4-2=2 but we have 3 ranks...
            // I need to remove the blocks that have no block from my send and myself if needed
            n_send_request_ = 0;
            for (rank_t ir = 0; ir < n_recver; ++ir) {
                const rank_t c_recver = first_recver + ir;
                if (c_recver == rank) {
                    continue;
                }
                n_send_request_ += (forest->global_first_quadrant[c_recver + 1] - forest->global_first_quadrant[c_recver]) > 0;
            }
            // if I am among the list, remove myself from the request point of view
            // n_send_request_ -= (first_recver <= rank && rank <= last_recver);
            m_verb("I will do %d send reqests to send blocks to the %d detected receivers", n_send_request_, n_recver);

            // allocate the requests
            for_send_request_  = reinterpret_cast<MPI_Request *>(m_calloc(n_send_request_ * sizeof(MPI_Request)));
            back_recv_request_ = reinterpret_cast<MPI_Request *>(m_calloc(n_send_request_ * sizeof(MPI_Request)));
            // allocate the cummulative list, to know which block I have to copy to the buffer
            // note: we cannot use a cummulative since we may not have continuous set (if some blocks are not moving)
            q_send_cum_block_   = reinterpret_cast<iblock_t *>(m_calloc(n_send_request_ * sizeof(iblock_t)));
            q_send_cum_request_ = reinterpret_cast<iblock_t *>(m_calloc((n_send_request_ + 1) * sizeof(iblock_t)));

            // for each receiver, allocate the send request
            rank_t   rcount  = 0;
            iblock_t qcount  = 0;
            iblock_t tqcount = 0;
            for (rank_t ir = 0; ir < n_recver; ++ir) {
                // get who is the receiver and skip if it's me
                const rank_t c_recver = first_recver + ir;
                // get howmany the receiver will receive (not his/her entire block numbers)
                const p4est_gloidx_t q_leftlimit  = m_max(forest->global_first_quadrant[c_recver], opart_begin);
                const p4est_gloidx_t q_rightlimit = m_min(forest->global_first_quadrant[c_recver + 1], opart_end);
                const iblock_t       n_q2send     = q_rightlimit - q_leftlimit;
                // if I want to allocate a send for myself, skip
                if (c_recver == rank || n_q2send == 0) {
                    tqcount += n_q2send;
                    continue;
                }
                // if we send the status, store the cumsum as well
                // if (destructive_) {
                //     m_assert(sizeof(bidx_t) == sizeof(int), "if not, we are wrong in the datatypes for mpi");
                //     send_status_count_[c_recver]   = n_q2send;
                //     send_status_cum_sum_[c_recver] = qcount;  // send_status_cum_sum_[c_recver] + n_q2send;
                // }
                // remember the begin and end point
                m_assert(rcount < n_send_request_, "scount = %d is too big compared to the number of request = %d", rcount, n_send_request_);
                q_send_cum_block_[rcount]   = tqcount;
                q_send_cum_request_[rcount] = qcount;

                // create the send request
                real_t *__restrict buf = send_buf_.ptr + qcount * partition_blocksize_;
                MPI_Send_init(buf, partition_blocksize_ * n_q2send, M_MPI_REAL, c_recver, c_recver, forest->mpicomm, &(for_send_request_[rcount]));
                MPI_Recv_init(buf, partition_blocksize_ * n_q2send, M_MPI_REAL, c_recver, c_recver, forest->mpicomm, &(back_recv_request_[rcount]));

                // update the counters
                qcount += n_q2send;
                tqcount += n_q2send;
                rcount += 1;
            }
            // store the last index
            q_send_cum_request_[n_send_request_] = qcount;
            m_assert(qcount == opart_n, "counters are not matching");

            m_verb("done with %d send requests", n_send_request_);
        } else {
            m_verb("No blocks to send");
        }

        //................................................
        // init the reception
        if (cpart_n > 0) {
            // store the new quadrant adress, to create a new block if needed
            new_blocks_                = reinterpret_cast<CartBlock **>(m_calloc(forest->local_num_quadrants * sizeof(CartBlock *)));
            p4est_gloidx_t rank_offset = forest->global_first_quadrant[rank];

            // got through each block and add it if needed
            iblock_t bcount = 0;
            for (p4est_topidx_t it = forest->first_local_tree; it <= forest->last_local_tree; it++) {
                p8est_tree_t *mytree = p8est_tree_array_index(forest->trees, it);
                for (iblock_t qid = 0; qid < mytree->quadrants.elem_count; qid++) {
                    qdrt_t * quad   = p8est_quadrant_array_index(&mytree->quadrants, qid);
                    iblock_t offset = mytree->quadrants_offset;
                    // get the global ID and allocate a new block only if the quad is new
                    p4est_gloidx_t global_id = rank_offset + offset + qid;
                    if (!(opart_begin <= global_id && global_id < opart_end)) {
                        // create a new block
                        real_t xyz[3];
                        p8est_qcoord_to_vertex(forest->connectivity, it, quad->x, quad->y, quad->z, xyz);
                        m_assert(quad->level >= 0, "the level=%d must be >=0", quad->level);
                        // real_t     len   = p4est_QuadLen(quad->level);
                        // CartBlock *block = AllocateGridBlock(grid->block_type(), len, xyz, quad->level);
                        AllocateBlockForP4est(grid->block_type(), xyz, quad);
                        CartBlock *block = p4est_GetBlock<CartBlock>(quad);
                        // add the fields to the block
                        block->AddFields(fields);
                        // store its access, replace the adress that was there but which is wrong now
                        // *(reinterpret_cast<GridBlock **>(quad->p.user_data)) = block;
                        // p4est_SetBlock(quad,block);
                        // store in the array
                        new_blocks_[offset + qid] = block;
                        // increment the counter
                        bcount++;
                    } else {
                        new_blocks_[offset + qid] = nullptr;  // *(reinterpret_cast<GridBlock **>(quad->p.user_data));
                    }
                }
            }
            m_assert(bcount == cpart_n, "the counters should match");

            // the receive buffer is used as a new data location for the blocks
            recv_buf_.Allocate(cpart_n * partition_blocksize_);
            // m_log("receiving buffer initialize of size %ld bytes (n_lda = %d)", cpart_n * PartCommSize(n_lda_) * sizeof(real_t), n_lda_);

            // if (destructive_) {
            //     // if destructive, send the status
            //     recv_status_buf_  = reinterpret_cast<short *>(m_calloc(cpart_n * sizeof(short)));
            // }

            // senders = from cpart_begin to cpart_end-1 in the new partition
            m_verb("looking for the senders of my new block: %ld -> %ld", cpart_begin, cpart_end);
            const rank_t first_sender = bsearch_comm(oldpart, cpart_begin, commsize, 0);
            const rank_t last_sender  = bsearch_comm(oldpart, cpart_end - 1, commsize, 0);
            const rank_t n_sender     = last_sender - first_sender + 1;  // note: add one because 4 - 2 = 2 but we have 3 ranks in this case
            m_verb("found sender: rank %d to %d", first_sender, last_sender);

            // we need to remove the ranks that have no blocks and myself if I am in the list
            n_recv_request_ = 0;
            for (rank_t ir = 0; ir < n_sender; ++ir) {
                const rank_t c_sender = first_sender + ir;
                if (c_sender == rank) {
                    continue;
                }
                n_recv_request_ += (oldpart[c_sender + 1] - oldpart[c_sender]) > 0;
            }
            // remove myself
            // n_recv_request_ -= (first_sender <= rank && rank <= last_sender);
            m_verb("I will do %d recv reqests to get blocks from %d senders", n_recv_request_, n_sender);

            // allocate the requests
            for_recv_request_  = reinterpret_cast<MPI_Request *>(m_calloc(n_recv_request_ * sizeof(MPI_Request)));
            back_send_request_ = reinterpret_cast<MPI_Request *>(m_calloc(n_recv_request_ * sizeof(MPI_Request)));
            // allocate the cummulative list, to know which block I have to copy to the buffer
            // note: we cannot use a cummulative since we may not have continuous set (if some blocks are not moving)
            q_recv_cum_block_   = reinterpret_cast<iblock_t *>(m_calloc(n_recv_request_ * sizeof(iblock_t)));
            q_recv_cum_request_ = reinterpret_cast<iblock_t *>(m_calloc((n_recv_request_ + 1) * sizeof(iblock_t)));

            // for each receiver, allocate the send request
            rank_t   scount  = 0;  // counter on the number of sender (rank)
            iblock_t qcount  = 0;  // counter on the blocks actually send
            iblock_t tqcount = 0;  // counter on all the blocks on this rank
            for (rank_t ir = 0; ir < n_sender; ir++) {
                // get who is the sender
                const rank_t c_sender = first_sender + ir;

                // get howmany the sender will send (not his/her entire block numbers)
                const p4est_gloidx_t q_leftlimit  = m_max(oldpart[c_sender], cpart_begin);
                const p4est_gloidx_t q_rightlimit = m_min(oldpart[c_sender + 1], cpart_end);
                const iblock_t       n_q2recv     = q_rightlimit - q_leftlimit;
                // if I am the sender, advance the tqcount counter and skip
                if (c_sender == rank || n_q2recv == 0) {
                    tqcount += n_q2recv;
                    continue;
                }
                // if we send the status, store the cumsum as well
                // if (destructive_) {
                //     m_assert(sizeof(bidx_t) == sizeof(int), "if not, we are wrong in the datatypes for mpi");
                //     recv_status_count_[c_sender]   = n_q2recv;
                //     recv_status_cum_sum_[c_sender] = qcount;
                // }

                // store the memory accesses
                m_assert(scount < n_recv_request_, "scount = %d is too big compared to the number of request = %d", scount, n_recv_request_);
                q_recv_cum_block_[scount]   = tqcount;
                q_recv_cum_request_[scount] = qcount;
                // create the send request
                real_t *__restrict buf = recv_buf_.ptr + qcount * partition_blocksize_;
                MPI_Recv_init(buf, partition_blocksize_ * n_q2recv, M_MPI_REAL, c_sender, rank, forest->mpicomm, &(for_recv_request_[scount]));
                MPI_Send_init(buf, partition_blocksize_ * n_q2recv, M_MPI_REAL, c_sender, rank, forest->mpicomm, &(back_send_request_[scount]));
                // update the counters
                qcount += n_q2recv;
                tqcount += n_q2recv;
                scount += 1;
            }
            q_recv_cum_request_[n_recv_request_] = qcount;
            m_assert(qcount == cpart_n, "counters are not matching");
            m_assert(qcount == cpart_n, "counters are not matching");
        } else {
            m_verb("No blocks to recv");
        }
    }
    m_free(oldpart);

    //................................................
    // reset the field list in the grid if we didn't conserve everything
    if (!destructive_) {
        grid->ResetFields(fields);
    }

    m_log("grid partitioned: moving %ld blocks (%f percent of the grid)", nqshipped, 100.0 * ((real_t)nqshipped) / ((real_t)forest->global_num_quadrants));

    //-------------------------------------------------------------------------
    m_end;
}

Partitioner::~Partitioner() {
    if (destructive_) {
        DeallocOldies_();
    }
    if (old_blocks_ != nullptr) {
        m_free(old_blocks_);
    }
    if (new_blocks_ != nullptr) {
        m_free(new_blocks_);
    }

    // if (destructive_) {
    //     m_free(send_status_count_);
    //     m_free(send_status_cum_sum_);
    //     m_free(recv_status_count_);
    //     m_free(recv_status_cum_sum_);
    //     m_free(send_status_buf_);
    //     m_free(recv_status_buf_);
    // }

    for (int i = 0; i < n_send_request_; i++) {
        MPI_Request_free(for_send_request_ + i);
        MPI_Request_free(back_recv_request_ + i);
    }
    for (int i = 0; i < n_recv_request_; i++) {
        MPI_Request_free(for_recv_request_ + i);
        MPI_Request_free(back_send_request_ + i);
    }

    if (n_send_request_ > 0) {
        send_buf_.Free();
    }

    if (n_recv_request_ > 0) {
        recv_buf_.Free();
    }

    m_free(for_send_request_);
    m_free(for_recv_request_);
    m_free(back_send_request_);
    m_free(back_recv_request_);
    m_free(q_send_cum_block_);
    m_free(q_send_cum_request_);
    m_free(q_recv_cum_block_);
    m_free(q_recv_cum_request_);
}

void Partitioner::SendRecv(map<string, Field *> *fields, const m_direction_t dir) {
    //--------------------------------------------------------------------------
    lid_t  n_recv_request;
    lid_t  n_send_request;
    lid_t *q_send_cum_request;
    lid_t *q_recv_cum_request;

    lid_t *q_send_cum_block;
    lid_t *q_recv_cum_block;

    real_t *send_buf;
    real_t *recv_buf;

    MPI_Request *recv_request;
    MPI_Request *send_request;
    CartBlock ** old_blocks;
    CartBlock ** new_blocks;

    if (dir == M_FORWARD) {
        n_recv_request     = n_recv_request_;
        n_send_request     = n_send_request_;
        q_send_cum_request = q_send_cum_request_;
        q_recv_cum_request = q_recv_cum_request_;
        q_send_cum_block   = q_send_cum_block_;
        q_recv_cum_block   = q_recv_cum_block_;
        send_buf           = send_buf_.ptr;
        recv_buf           = recv_buf_.ptr;
        recv_request       = for_recv_request_;
        send_request       = for_send_request_;
        new_blocks         = new_blocks_;
        old_blocks         = old_blocks_;
    } else {
        n_recv_request     = n_send_request_;
        n_send_request     = n_recv_request_;
        q_send_cum_request = q_recv_cum_request_;
        q_recv_cum_request = q_send_cum_request_;
        q_send_cum_block   = q_recv_cum_block_;
        q_recv_cum_block   = q_send_cum_block_;
        send_buf           = recv_buf_.ptr;
        recv_buf           = send_buf_.ptr;
        recv_request       = back_recv_request_;
        send_request       = back_send_request_;
        old_blocks         = new_blocks_;
        new_blocks         = old_blocks_;

        m_assert(!destructive_, "we cannot send backward and be destructive, otherwise the status are not sent");
    }

    //..........................................................................
    // performs the send of the request is
    auto send_my_request = [=](const lid_t is) {
        lid_t n_q2send = q_send_cum_request[is + 1] - q_send_cum_request[is];
        for (lid_t iq = 0; iq < n_q2send; iq++) {
            // get the block to send
            CartBlock *block = old_blocks[q_send_cum_block[is] + iq];
            real_t *   buf   = send_buf + (q_send_cum_request[is] + iq) * partition_blocksize_;
            m_assert(block != nullptr, "this block shouldn't be nullptr here");
            m_assert(partition_blocksize_ == PartitionerBlockSize(block, n_lda_), "the two sizes must match: %ld vs %ld", partition_blocksize_, PartitionerBlockSize(block, n_lda_));

            // the first numbers are the status
            block->PartitionDataPack(buf);
            buf += block->PartitionDataOffset();

            // handle the field
            bidx_t       idacount   = 0;
            const size_t block_size = block->BlockLayout().n_elem;
            for (const auto iter : *fields) {
                const Field *fid  = iter.second;
                real_t *     lbuf = buf + block_size * idacount;
                // copy the whole memory at once on the dim
                if (!fid->is_expr()) {
                    memcpy(lbuf, block->RawPointer(fid, 0), block_size * fid->lda() * sizeof(real_t));
                    // update the counters
                    idacount += fid->lda();
                }

                m_assert(fid->lda() <= n_lda_, "unable to transfert so much data with the allocated buffers");
            }
        }
        // start the send
        MPI_Start(send_request + is);
        m_verb("I have started the send of request %d", is);
    };

    // performs the recv of the request idx
    auto recv_my_request = [=](const lid_t ir) {
        lid_t n_q2recv = q_recv_cum_request[ir + 1] - q_recv_cum_request[ir];
        for (lid_t iq = 0; iq < n_q2recv; iq++) {
            // get the block
            CartBlock *block = new_blocks[q_recv_cum_block[ir] + iq];
            real_t *   buf   = recv_buf + (q_recv_cum_request[ir] + iq) * partition_blocksize_;
            m_assert(block != nullptr, "this block shouldn't be nullptr here");
            m_assert(partition_blocksize_ == PartitionerBlockSize(block, n_lda_), "the two sizes must match: %ld vs %ld", partition_blocksize_, PartitionerBlockSize(block, n_lda_));

            // the first numbers are the status
            block->PartitionDataUnPack(buf);
            buf += block->PartitionDataOffset();

            // handle the field
            bidx_t       idacount   = 0;
            const size_t block_size = block->BlockLayout().n_elem;
            for (auto iter : *fields) {
                const Field *fid = iter.second;
                // security check
                m_assert(fid->lda() <= n_lda_, "unable to transfert so much data with the allocated buffers");
                const real_t *lbuf = buf + block_size * idacount;
                // copy the whole memory at once on the dim
                if (!fid->is_expr()) {
                    memcpy(block->RawPointer(fid, 0), lbuf, block_size * fid->lda() * sizeof(real_t));
                    // update the counters
                    idacount += fid->lda();
                }
            }
        }
        m_verb("I have finished the recv of request %d", ir);
    };

    //..........................................................................
    // we do the send while waiting for a recv request to complete
    // start all the reception
    if (n_recv_request > 0) {
        MPI_Startall(n_recv_request, recv_request);
    }

    // allocate the worst case array
    int * current_recv_id = reinterpret_cast<int *>(m_calloc(n_recv_request * sizeof(int)));
    lid_t n_recv_done     = 0;
    lid_t n_send_done     = 0;
    while ((n_recv_done < n_recv_request) || (n_send_done < n_send_request)) {
        // m_log("already done %d recv and %d sends", n_recv_done, n_send_done);
        // if we have some recv left, test them
        if (n_recv_done < n_recv_request) {
            int n_recv_ready = 0;
            MPI_Testsome(n_recv_request, recv_request, &n_recv_ready, current_recv_id, MPI_STATUSES_IGNORE);
            m_assert(MPI_UNDEFINED < 0, "the MPI undefined value must be <0");
            for (int ir = 0; ir < n_recv_ready; ++ir) {
                recv_my_request(current_recv_id[ir]);
                n_recv_done += 1;
            }
        }
        // once done, send 1 request
        if (n_send_done < n_send_request) {
            send_my_request(n_send_done);
            n_send_done += 1;
        }
    }
    m_assert((n_send_done == n_send_request) && (n_recv_done == n_recv_request), "the number of requests done must = the number of total requests: %d vs %d and %d vs %d", n_send_done, n_send_request, n_recv_done, n_recv_request);
    m_free(current_recv_id);

    // wait for all the send to be ok (should be trivial)
    if (n_send_request > 0) {
        MPI_Waitall(n_send_request, send_request, MPI_STATUSES_IGNORE);
    }

#ifndef NDEBUG
    // no need to Waitall on the recv requests as they have been processed
    // if in debug stil do a test on all the requests
    for (lid_t is = 0; is < n_send_request; ++is) {
        // get the status of the request and check it's over
        int        is_over = 0;
        MPI_Status status;
        MPI_Request_get_status(send_request[is], &is_over, &status);
        m_assert(is_over, "The send request num %d with source %d is not finished from ", is, status.MPI_SOURCE);
    }
    for (lid_t ir = 0; ir < n_recv_request; ++ir) {
        // get the status of the request and check it's over
        int        is_over = 0;
        MPI_Status status;
        MPI_Request_get_status(recv_request[ir], &is_over, &status);
        m_assert(is_over, "The send request num %d with source %d is not finished from ", ir, status.MPI_SOURCE);
    }
#endif
}

/**
 * @brief starts the send communication from the old paritioning to the new one: copy the GridBlock to the buffer and start the send and recv
 * 
 * @warning the ghost region is sent with the data for two reasons: (1) it avoids weird buffer as they are not continuous or MPI_Datatypes and (2)
 * we need the ghost values to solve the resolution jump, which is happening after the partitioning
 * 
 * @param fields the field(s) that will be transfered (only one field is accepted if the non-destructive mode is enabled)
 * @param dir the direction, forward or backward
 */
void Partitioner::Start(map<string, Field *> *fields, const m_direction_t dir) {
    m_begin;
    m_assert(dir == M_FORWARD || dir == M_BACKWARD, "the dir can be forward or backward");
    m_assert(!(dir == M_BACKWARD && destructive_), "unable to perform backward on destructive partitioner");
    m_assert(!(fields->size() > 1 && !destructive_), "the partitioner has been allocated in destructive mode, unable to transfert mode than 1 field");
    //-------------------------------------------------------------------------
    rank_t    n_recv_request;
    rank_t    n_send_request;
    iblock_t *q_send_cum_request;
    iblock_t *q_send_cum_block;
    real_t *  send_buf;

    MPI_Request *recv_request;
    MPI_Request *send_request;
    CartBlock ** old_blocks;

    if (dir == M_FORWARD) {
        n_recv_request     = n_recv_request_;
        n_send_request     = n_send_request_;
        recv_request       = for_recv_request_;
        send_request       = for_send_request_;
        send_buf           = send_buf_.ptr;
        old_blocks         = old_blocks_;
        q_send_cum_block   = q_send_cum_block_;
        q_send_cum_request = q_send_cum_request_;
    } else {
        n_recv_request     = n_send_request_;
        n_send_request     = n_recv_request_;
        recv_request       = back_recv_request_;
        send_request       = back_send_request_;
        send_buf           = recv_buf_.ptr;
        old_blocks         = new_blocks_;
        q_send_cum_block   = q_recv_cum_block_;
        q_send_cum_request = q_recv_cum_request_;

        m_assert(!destructive_, "we cannot send backward and be destructive, otherwise the status are not sent");
    }

    // start the reception of the memory
    if (n_recv_request > 0) {
        MPI_Startall(n_recv_request, recv_request);
    }
    // copy all the data of the old blocks into the send_buffer and send when ready
    for (rank_t is = 0; is < n_send_request; ++is) {
        iblock_t n_q2send = q_send_cum_request[is + 1] - q_send_cum_request[is];
        for (rank_t iq = 0; iq < n_q2send; ++iq) {
            const CartBlock *block = old_blocks[q_send_cum_block[is] + iq];
            real_t *         buf   = send_buf + (q_send_cum_request[is] + iq) * partition_blocksize_;
            m_assert(block != nullptr, "this block shouldn't be nullptr here");
            m_assert(partition_blocksize_ == PartitionerBlockSize(block, n_lda_), "the two sizes must match: %ld vs %ld", partition_blocksize_, PartitionerBlockSize(block, n_lda_));

            // the first numbers are the status
            block->PartitionDataPack(buf);
            buf += block->PartitionDataOffset();

            // get the data now
            const size_t block_size = block->BlockLayout().n_elem;
            lda_t        idacount   = 0;
            for (const auto iter : *fields) {
                const Field *fid = iter.second;
                // security check
                m_assert(fid->lda() <= n_lda_, "unable to transfert so much data with the allocated buffers");
                // real_t* data = ;
                real_t *lbuf = buf + block_size * idacount;
                // copy the whole memory at once on the dim
                if (!fid->is_expr()) {
                    memcpy(lbuf, block->RawPointer(fid, 0), block_size * fid->lda() * sizeof(real_t));
                    // update the counters
                    idacount += fid->lda();
                }
            }
        }
        // start the send
        MPI_Start(send_request + is);
    }
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief terminates the partitioning: end the reception, the send and copy the buffer content to the new blocks
 * 
 * @param fields the fields that will be transfered (may not correspond to the initialization one but HAS to match the one used in the Start())
 * @param dir the direction, forward or backward
 * @param do_copy do we need to copy the non-traveling blocks or not
 */
void Partitioner::End(map<string, Field *> *fields, const m_direction_t dir) {
    m_begin;
    m_assert(dir == M_FORWARD || dir == M_BACKWARD, "the dir can be forward or backward");
    m_assert(!(dir == M_BACKWARD && destructive_), "unable to perform backward on destructive partitioner");
    m_assert(!(fields->size() > 1 && !destructive_), "the partitioner has been allocated in destructive mode, unable to transfert mode than 1 field");
    //-------------------------------------------------------------------------
    rank_t    n_recv_request;
    rank_t    n_send_request;
    iblock_t *q_recv_cum_request;
    iblock_t *q_recv_cum_block;
    real_p    recv_buf;

    MPI_Request *recv_request;
    MPI_Request *send_request;
    CartBlock ** new_blocks;

    if (dir == M_FORWARD) {
        n_recv_request     = n_recv_request_;
        n_send_request     = n_send_request_;
        q_recv_cum_request = q_recv_cum_request_;
        q_recv_cum_block   = q_recv_cum_block_;
        recv_buf           = recv_buf_.ptr;
        recv_request       = for_recv_request_;
        send_request       = for_send_request_;
        new_blocks         = new_blocks_;
    } else {
        n_recv_request     = n_send_request_;
        n_send_request     = n_recv_request_;
        q_recv_cum_request = q_send_cum_request_;
        q_recv_cum_block   = q_send_cum_block_;
        recv_buf           = send_buf_.ptr;
        recv_request       = back_recv_request_;
        send_request       = back_send_request_;
        new_blocks         = old_blocks_;
        m_assert(!destructive_, "we cannot send backward and be destructive, otherwise the status are not sent");
    }

    // handle the moving ones
    for (rank_t ir = 0; ir < n_recv_request; ++ir) {
        rank_t idx;
        MPI_Waitany(n_recv_request, recv_request, &idx, MPI_STATUS_IGNORE);
        iblock_t n_q2recv = q_recv_cum_request[idx + 1] - q_recv_cum_request[idx];

        for (iblock_t iq = 0; iq < n_q2recv; ++iq) {
            CartBlock *block = new_blocks[q_recv_cum_block[idx] + iq];
            real_t *   buf   = recv_buf + (q_recv_cum_request[idx] + iq) * partition_blocksize_;
            m_assert(block != nullptr, "this block shouldn't be nullptr here");
            m_assert(partition_blocksize_ == PartitionerBlockSize(block, n_lda_), "the two sizes must match: %ld vs %ld", partition_blocksize_, PartitionerBlockSize(block, n_lda_));

            // the first numbers are the status
            block->PartitionDataUnPack(buf);
            buf += block->PartitionDataOffset();

            bidx_t       idacount   = 0;
            const size_t block_size = block->BlockLayout().n_elem;
            for (auto iter : *fields) {
                const Field *fid = iter.second;
                // security check
                m_assert(fid->lda() <= n_lda_, "unable to transfert so much data with the allocated buffers");
                const real_t *lbuf = buf + block_size * idacount;
                // copy the whole memory at once on the dim
                if (!fid->is_expr()) {
                    memcpy(block->RawPointer(fid, 0), lbuf, block_size * fid->lda() * sizeof(real_t));
                    // update the counters
                    idacount += fid->lda();
                }
            }
        }
    }
    // receive the new memory in the new blocks
    if (n_send_request > 0) {
        MPI_Waitall(n_send_request, send_request, MPI_STATUSES_IGNORE);
    }
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief deallocates the old GridBlock if they are not needed anymore
 * 
 */
void Partitioner::DeallocOldies_() {
    m_begin;
    //-------------------------------------------------------------------------
    for (int is = 0; is < n_send_request_; is++) {
        // delete all the blocks that have been transfered, not the remaining ones..
        iblock_t n_q2send = q_send_cum_request_[is + 1] - q_send_cum_request_[is];
        for (iblock_t iq = q_send_cum_block_[is]; iq < (q_send_cum_block_[is] + n_q2send); iq++) {
            if (old_blocks_[iq] != nullptr) {
                delete (old_blocks_[iq]);
            }
        }
    }
    //-------------------------------------------------------------------------
    m_end;
}
