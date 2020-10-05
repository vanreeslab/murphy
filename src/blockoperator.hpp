#ifndef SRC_BLOCKOPERATOR_HPP
#define SRC_BLOCKOPERATOR_HPP

#include "gridblock.hpp"
#include "murphy.hpp"
#include "interpolator.hpp"

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

   public:
    /**
    * @brief Construct a new Block Operator object
    * 
    * @param interp if the interp is nullptr, no ghost point is filled, if not, we fill it
    */
    explicit BlockOperator(const Interpolator* interp) {
        m_begin;
        //-------------------------------------------------------------------------
        do_ghost_ = (interp != nullptr);
        start_    = (interp == nullptr) ? 0 : (-interp->nghost_front());
        end_      = (interp == nullptr) ? M_N : (M_N + interp->nghost_back());
        //-------------------------------------------------------------------------
        m_end;
    }

    /**
     * @brief return true if the BlockOperator fills the ghosts
     */
    inline bool do_ghost() const { return do_ghost_; }
};

#endif