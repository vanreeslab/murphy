#ifndef SRC_BLOCKOPERATOR_HPP
#define SRC_BLOCKOPERATOR_HPP

#include "core/memspan.hpp"
#include "core/macros.hpp"
#include "core/types.hpp"
#include "tools/prof.hpp"

/**
 * @brief represents the simplest operator on a block that operates from [start_,start_,start_] to [end_,end_,end_]
 * 
 * No structure is imposed for the operator, the call or anything else.
 * We just take into account the simplest redundant operation: the start and end index and the storage of the profiler
 * 
 */
class BlockOperator {
   protected:
    const MemSpan span_;            //!< Span for the loop
    
    bidx_t       ghost_len_need_[2] = {0, 0};   //!< the number of ghost points needed by the operator
    bidx_t       ghost_len_res_[2]  = {0, 0};   //!< the number of ghost points actually visited by the operator and ok in the result
    Prof*        prof_              = nullptr;  //!< profiler

   public:
    explicit BlockOperator() = delete;
    explicit BlockOperator(const bidx_t* ghost_len_in) noexcept : BlockOperator(ghost_len_in, ghost_len_in){};

    explicit BlockOperator(const bidx_t* ghost_len_need, const bidx_t* ghost_len_res) noexcept
        : span_(((ghost_len_res == nullptr) ? 0 : (-ghost_len_res[0])),
          ((ghost_len_res == nullptr) ? M_N : (M_N + ghost_len_res[1]))) {
        m_begin;
        //-------------------------------------------------------------------------
        ghost_len_need_[0] = (ghost_len_need == nullptr) ? 0 : ghost_len_need[0];
        ghost_len_need_[1] = (ghost_len_need == nullptr) ? 0 : ghost_len_need[1];
        ghost_len_res_[0]  = (ghost_len_res == nullptr) ? 0 : ghost_len_res[0];
        ghost_len_res_[1]  = (ghost_len_res == nullptr) ? 0 : ghost_len_res[1];
        //-------------------------------------------------------------------------
        m_end;
    };

    [[nodiscard]] inline bool DoGhost() const { return (ghost_len_res_[0] + ghost_len_res_[1]) > 0; }

    inline virtual void Profile(Prof* profiler) { prof_ = profiler; }

    inline void GhostLengthNeed(bidx_t* ghost_len) const {
        ghost_len[0] = ghost_len_need_[0];
        ghost_len[1] = ghost_len_need_[1];
    }
    inline void GhostLenghtResult(bidx_t* ghost_len) const {
        ghost_len[0] = ghost_len_res_[0];
        ghost_len[1] = ghost_len_res_[1];
    }
};

#endif