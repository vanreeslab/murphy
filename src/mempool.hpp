#ifndef MEMPOOL_HPP_
#define MEMPOOL_HPP_

#include <list>

#include "murphy.hpp"

using std::list;

/**
 * @brief a thread safe pooling for temporary/working memory
 * 
 */
class MemPool {
   private:
    size_t         mem_size_ = 0;  //!< the size of a memory block allocated
    list<real_t*>* mem_pool_;      //!< the list of the temporary memory allocated

   public:
    explicit MemPool();
    explicit MemPool(size_t mem_size);
    ~MemPool();

    real_t* LockMemory(sid_t ithread);
    void    FreeMemory(sid_t ithread, real_t*);
};

#endif  // MEMPOOL_HPP_