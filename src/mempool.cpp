#include "mempool.hpp"

#include <limits>
#include <string>

using std::numeric_limits;
using std::string;

MemPool::MemPool():MemPool(m_blockmemsize(1)) {}

MemPool::MemPool(size_t mem_size) {
    m_begin;
    //-------------------------------------------------------------------------
    // remember the memsize
    mem_size_ = mem_size;
    // get the number of max threads and check that not too high
    m_assert(omp_get_max_threads() < numeric_limits<sid_t>::max(), "the number of threads is too high");
    sid_t nthread = omp_get_max_threads();

    // allocate the pool, one per thread
    n_locked_ = reinterpret_cast<lid_t*>(m_calloc(sizeof(lid_t) * nthread));
    mem_pool_ = reinterpret_cast<list<real_t*>**>(m_calloc(sizeof(list<real_t*>*) * nthread));
    for (sid_t it = 0; it < nthread; it++) {
        n_locked_[it] = 0;
        mem_pool_[it] = new list<real_t*>();
    }

    //-------------------------------------------------------------------------
    m_end;
}

MemPool::~MemPool() {
    m_begin;
    //-------------------------------------------------------------------------
    sid_t nthread = omp_get_max_threads();
    for (sid_t it = 0; it < nthread; it++) {
        m_assert(n_locked_[it] == 0, "the number of locked array has to be 0 before the exit");
        // deallocate the allocated blocks
        for(auto iter = mem_pool_[it]->begin(); iter!= mem_pool_[it]->end(); iter++){
            m_free(*iter);
        }
        // clear the list
        mem_pool_[it]->clear();
    }
    m_free(n_locked_);
    m_free(mem_pool_);
    //-------------------------------------------------------------------------
    m_end;
}

real_t* MemPool::LockMemory(sid_t ithread) {
    //-------------------------------------------------------------------------
    real_t* ptr;
    if (mem_pool_[ithread]->empty()) {
        ptr = reinterpret_cast<real_t*>(m_calloc(sizeof(real_t) * mem_size_));
        n_locked_++;
        m_log("allocating new memory of size %d",mem_size_);
    } else {
        ptr = mem_pool_[ithread]->back();
        memset(ptr, 0, sizeof(real_t) * mem_size_);
        mem_pool_[ithread]->pop_back();
        n_locked_++;
    }
    return ptr;
    //-------------------------------------------------------------------------
}

void MemPool::FreeMemory(sid_t ithread, real_t* ptr) {
    //-------------------------------------------------------------------------
    mem_pool_[ithread]->push_back(ptr);
    n_locked_--;
    //-------------------------------------------------------------------------
}