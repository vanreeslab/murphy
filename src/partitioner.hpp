#ifndef SRC_PARTITIONER_HPP_
#define SRC_PARTITIONER_HPP_

#include "murphy.hpp"
#include "forestgrid.hpp"
#include "gridblock.hpp"

class Partitioner {
   protected:
   lid_t n_lda_  = 0;
    lid_t        n_send_request_ = 0;
    lid_t        n_recv_request_ = 0;
    MPI_Request* send_request_   = nullptr;
    MPI_Request* recv_request_   = nullptr;

    real_p send_buf_ = nullptr;
    real_p recv_buf_ = nullptr;

    
    lid_t* q_send_start_ = nullptr;
    lid_t* q_send_end_ = nullptr;
    GridBlock** old_blocks_;

    //note: since every block owns his own memory, I have to copy it to the buffer
    lid_t* q_recv_start_ = nullptr;
    lid_t* q_recv_end_ = nullptr;
    GridBlock** new_blocks_;

   public:
    Partitioner(map<string, Field*>* fields, ForestGrid *grid) ;
    ~Partitioner();

    void Start(map<string, Field*>* fields);
    void End(map<string, Field*>* fields);

    protected:
    void DeallocOldies_();

};

#endif  // SRC_PARTITIONER_HPP_