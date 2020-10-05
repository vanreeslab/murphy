#ifndef SRC_BLOCKOPERATOR_HPP
#define SRC_BLOCKOPERATOR_HPP

#include "gridblock.hpp"
#include "murphy.hpp"
#include "interpolator.hpp"
#include "prof.hpp"

/**
 * @brief represents the simplest operator on a GridBlock that operates from [start_,start_,start_] to [end_,end_,end_] on a GridBlock
 * 
 * No structure is imposed for the operator, the call or anything else.
 * We just take into account the simplest redundant operation: the start and end index
 * 
 * @tparam T the template to handle any combination of the arguments, must be at the end of the call
 */
class BlockOperator {
   protected:
    bool  do_ghost_ = false;
    lid_t start_    = 0;
    lid_t end_      = M_N;

    Prof* prof_ = nullptr;

   public:
   /**
    * @brief cfr BlockOperator::BlockOperator(const Interpolator* interp, Prof* profiler)
    */
    // explicit BlockOperator() : BlockOperator(nullptr, nullptr) {}
    // explicit BlockOperator(const Interpolator* interp) : BlockOperator(interp, nullptr) {}
    // explicit BlockOperator(Prof* profiler) : BlockOperator(nullptr, profiler) {}

    /**
     * @brief Construct a new Block Operator object to operate based on the interpolation ghost points. The execution is timed using the profiler
     * 
     * @param interp the interpolation driving the number of ghost point to operate on. If nullptr, we operates from 0 to M_N
     * @param profiler the profier
     */
    explicit BlockOperator(const Interpolator* interp, Prof* profiler) {
        m_begin;
        //-------------------------------------------------------------------------
        do_ghost_ = (interp != nullptr);
        start_    = (interp == nullptr) ? 0 : (-interp->nghost_front());
        end_      = (interp == nullptr) ? M_N : (M_N + interp->nghost_back());
        prof_     = (profiler == nullptr) ? nullptr : profiler;
        //-------------------------------------------------------------------------
        m_end;
    }

    /**
     * @brief return true if the BlockOperator fills the ghosts
     */
    inline bool do_ghost() const { return do_ghost_; }
};

#endif