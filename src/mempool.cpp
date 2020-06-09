#include "mempool.hpp"

#include <limits>
#include <string>

using std::numeric_limits;
using std::string;

MemPool::MemPool() {
    MemPool(m_blockmemsize(1));
}

MemPool::MemPool(size_t mem_size) {
    m_begin;
    //-------------------------------------------------------------------------
    // remember the memsize
    mem_size_ = mem_size;
    // get the number of max threads and check that not too high
    m_assert(omp_get_max_threads() < numeric_limits<sid_t>::max(), "the number of threads is too high");
    sid_t nthread = omp_get_max_threads();

    // allocate the pool, one per thread
    mem_pool_ = reinterpret_cast<list<real_t*>*>(m_calloc(sizeof(list<real_t*>) * nthread));
    for (sid_t it = 0; it < nthread; it++) {
        mem_pool_[it] = list<real_t*>();
    }

    //-------------------------------------------------------------------------
    m_end;
}

real_t* MemPool::LockMemory(sid_t ithread) {
    //-------------------------------------------------------------------------
    real_t* ptr;
    if (mem_pool_[ithread].empty()) {
        ptr = reinterpret_cast<real_t*>(m_calloc(sizeof(real_t) * mem_size_));
    } else {
        ptr = mem_pool_[ithread].back();
        memset(ptr, 0, sizeof(real_t) * mem_size_);
        mem_pool_[ithread].pop_back();
    }
    return ptr;
    //-------------------------------------------------------------------------
}

void MemPool::FreeMemory(sid_t ithread, real_t* ptr) {
    //-------------------------------------------------------------------------
    mem_pool_[ithread].push_back(ptr);
    //-------------------------------------------------------------------------
}