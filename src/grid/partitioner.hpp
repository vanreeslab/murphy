#ifndef SRC_PARTITIONER_HPP_
#define SRC_PARTITIONER_HPP_

#include <map>
#include <string>

#include "core/macros.hpp"
#include "core/pointers.hpp"
#include "core/types.hpp"
#include "grid/grid.hpp"
#include "grid/gridblock.hpp"

/**
 * @brief Paritions a @ref Grid by distributing it on the available processors and send/recv the associated GridBlock*
 * 
 * The partitioner has to 'modes' dependings on the value of @ref destructive
 * If destructive is true: we do NOT care about the previous grid and it must be destroyed.
 *      - when destroying the Partitioner, it will deallocate the old blocks
 *      - the fields allocated must be all the fields of the Grid
 *      - the communication is done all the fields at once
 *      - the status_level_ of the block is transfered to the new block
 * If destructive is false: we DO care about the previous grid and it cannot be destroyed.
 *      - when destroying the partitioner, the old blocks are kept untouched
 *      - we can restrict the fields of the Grid to not allocate all of them
 *      - the communication is done one field at a time and hence requires less memory
 *      - the status_level_ of the block is NOT transfered to the new block
 */
class Partitioner {
   protected:
    bool  destructive_ = false;  //!< false if the oldblocks belongs to an existing grid
    lda_t n_lda_       = 0;      //!< the total number of dimension accross all the fields

    lid_t n_send_request_ = 0;  //!< the number of send requests
    lid_t n_recv_request_ = 0;  //!< the number of receive requests

    MPI_Request* for_send_request_  = nullptr;  //!< forward partition - the send requests
    MPI_Request* for_recv_request_  = nullptr;  //!< forward partition - the receive requests
    MPI_Request* back_send_request_ = nullptr;  //!< backward partition - the send requests
    MPI_Request* back_recv_request_ = nullptr;  //!< backward partition - the receive requests

    mem_ptr send_buf_ = nullptr;  //<! the send buffer, since the memory is not continuous accross the blocks
    mem_ptr recv_buf_ = nullptr;  //<! the receive buffer, sicne the memory is not continuous accross the blocks

    int*   send_status_count_   = nullptr;  //<! count to send the status of the old blocks
    int*   send_status_cum_sum_ = nullptr;  //<! cummulative sum to send the status of the old blocks
    int*   recv_status_count_   = nullptr;  //<! count to recv the status of the new blocks
    int*   recv_status_cum_sum_ = nullptr;  //<! cummulative sum to recv the status of the new blocks
    short* send_status_buf_     = nullptr;  //<! the send buffer for the status, since the memory is not continuous accross the blocks
    short* recv_status_buf_     = nullptr;  //<! the receive buffer for the status, sicne the memory is not continuous accross the blocks

    iblock_t*   q_send_cum_block_   = nullptr;  //!< remember at which local block we started for the request[i]
    iblock_t*   q_send_cum_request_ = nullptr;  //!< cummulative count on the moving blocks for the request[i]
    GridBlock** old_blocks_         = nullptr;  //!< array with the addresses to each old block

    // note: since every block owns his own memory, I have to copy it to the buffer for the receive as well
    iblock_t*   q_recv_cum_block_   = nullptr;  //!< remember at which local block we started for the request[i]
    iblock_t*   q_recv_cum_request_ = nullptr;  //!< cummulative count on the moving blocks for the request[i]
    GridBlock** new_blocks_         = nullptr;

   public:
    Partitioner(std::map<std::string, m_ptr<Field> >* fields, Grid* grid, bool destructive);
    ~Partitioner();

    void Start(std::map<std::string, m_ptr<Field> >* fields, const m_direction_t dir);
    void End(std::map<std::string, m_ptr<Field> >* fields, const m_direction_t dir);

   protected:
    void DeallocOldies_();
};

#endif  // SRC_PARTITIONER_HPP_
