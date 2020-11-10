#ifndef SRC_BLOCKOPERATOR_HPP
#define SRC_BLOCKOPERATOR_HPP

#include "gridblock.hpp"
#include "murphy.hpp"
#include "Wavelet.hpp"
#include "prof.hpp"

/**
 * @brief represents the simplest operator on a GridBlock that operates from [start_,start_,start_] to [end_,end_,end_] on a GridBlock
 * 
 * No structure is imposed for the operator, the call or anything else.
 * We just take into account the simplest redundant operation: the start and end index and the storage of the profiler
 * 
 */
class BlockOperator {
   protected:
    bool  do_ghost_ = false;
    lid_t start_    = 0;
    lid_t end_      = M_N;
    Prof* prof_     = nullptr;

   public:
    /**
     * @brief Construct a new Block Operator object to operate based on the interpolation ghost points. The execution is timed using the profiler
     * 
     * @param interp the interpolation driving the number of ghost point to operate on. If nullptr, we operates from 0 to M_N
     * @param profiler the profier
     */
    explicit BlockOperator(const Wavelet* interp) {
        m_begin;
        //-------------------------------------------------------------------------
        do_ghost_ = (interp != nullptr);
        start_    = (interp == nullptr) ? 0 : (-interp->nghost_front());
        end_      = (interp == nullptr) ? M_N : (M_N + interp->nghost_back());
        //-------------------------------------------------------------------------
        m_end;
    }

    void Profile(Prof* profiler) { prof_ = profiler; }

    /**
     * @brief return true if the BlockOperator fills the ghosts
     */
    inline bool do_ghost() const { return do_ghost_; }
};

#endif