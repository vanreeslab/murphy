#ifndef SRC_BLOCKOPERATOR_HPP
#define SRC_BLOCKOPERATOR_HPP

#include "core/macros.hpp"
#include "core/types.hpp"
#include "grid/gridblock.hpp"
#include "prof.hpp"
#include "wavelet/wavelet.hpp"

/**
 * @brief represents the simplest operator on a GridBlock that operates from [start_,start_,start_] to [end_,end_,end_] on a GridBlock
 * 
 * No structure is imposed for the operator, the call or anything else.
 * We just take into account the simplest redundant operation: the start and end index and the storage of the profiler
 * 
 */
class BlockOperator {
   protected:
    bool        do_ghost_ = false;    //!< indicate if the block operator fills the ghost region or not.
    lid_t       start_    = 0;        //!< starting index for the loop
    lid_t       end_      = M_N;      //!< ending index for the loop
    m_ptr<Prof> prof_     = nullptr;  //!< profiler

   public:
    /**
     * @brief Construct a new Block Operator object to operate based on the interpolation ghost points. The execution is timed using the profiler
     * 
     * @param interp the interpolation driving the number of ghost point to operate on. If nullptr, we operates from 0 to M_N
     * @param profiler the profier
     */
    explicit BlockOperator(m_ptr<const Wavelet> interp) {
        m_begin;
        //-------------------------------------------------------------------------
        do_ghost_ = !(interp.IsEmpty());
        start_    = (interp.IsEmpty()) ? 0 : (-interp()->nghost_front());
        end_      = (interp.IsEmpty()) ? M_N : (M_N + interp()->nghost_back());
        //-------------------------------------------------------------------------
        m_end;
    }

    void Profile(m_ptr<Prof> profiler) { prof_ = profiler; }

    /**
     * @brief return true if the BlockOperator fills the ghosts
     */
    inline bool do_ghost() const { return do_ghost_; }
};

#endif