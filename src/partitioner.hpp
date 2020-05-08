#ifndef SRC_PARTITIONER_HPP_
#define SRC_PARTITIONER_HPP_

#include <string>

#include "forestgrid.hpp"
#include "gridblock.hpp"
#include "murphy.hpp"

/**
 * @brief Paritions a ForestGrid by distributing it on the available processors and send/recv the associated GridBlock*
 * 
 */
class Partitioner {
   protected:
    bool         dangling_oldies_ = false;    //!< false if the oldblocks belongs to an existing grid
    lid_t        n_lda_           = 0;        //!< the total number of dimension accross all the fields
    lid_t        n_send_request_  = 0;        //!< the number of send requests
    lid_t        n_recv_request_  = 0;        //!< the number of receive requests
    MPI_Request* send_request_    = nullptr;  //!< the send requests
    MPI_Request* recv_request_    = nullptr;  //!< the receive requests

    real_p send_buf_ = nullptr;  //<! the send buffer, since the memory is not continuous accross the blocks
    real_p recv_buf_ = nullptr;  //<! the receive buffer, sicne the memory is not continuous accross the blocks

    lid_t*      q_send_cum_block_ = nullptr;  //!< remember at which local block we started for the request[i]
    lid_t*      q_send_cum_request_ = nullptr;  //!< cummulative count on the moving blocks for the request[i]
    GridBlock** old_blocks_         = nullptr;

    // note: since every block owns his own memory, I have to copy it to the buffer for the receive as well
    lid_t*      q_recv_cum_block_ = nullptr;  //!< remember at which local block we started for the request[i]
    lid_t*      q_recv_cum_request_ = nullptr;  //!< cummulative count on the moving blocks for the request[i]
    GridBlock** new_blocks_         = nullptr;

   public:
    Partitioner(map<string, Field*>* fields, ForestGrid* grid,bool dangling);
    ~Partitioner();

    void Start(map<string, Field*>* fields);
    void End(map<string, Field*>* fields);

   protected:
    void DeallocOldies_();
};

#endif  // SRC_PARTITIONER_HPP_
