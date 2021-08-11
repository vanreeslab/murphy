#ifndef SRC_CORE_MEMSPAN_HPP_
#define SRC_CORE_MEMSPAN_HPP_

#include "core/types.hpp"
#include "core/macros.hpp"

typedef enum m_layout_t {
    M_LAYOUT_BLOCK = 3
} m_layout_type_t;

struct MemLayout {
    size_t gs;         //!< number of ghost points in front of the block
    size_t shift;      //!< the shift to apply to the fastest rotating dimension to have the 0 aligned
    size_t stride[2];  //!< the strides of the data: fastest rotating stride in [0], user stride in [1]

    size_t n_elem;  //!< number of elements in the memory layout (= minimal allocation size)

    explicit MemLayout() = delete;
    explicit MemLayout(const m_layout_t layout, const bidx_t n_gs_front, const bidx_t n_block, const bidx_t n_gs_back = -1) noexcept;

    __attribute__((always_inline)) inline bidx_t offset(const bidx_t i0, const bidx_t i1, const bidx_t i2) const noexcept {
        return shift + (gs+i0) + stride[0] * ((gs+i1) + stride[1] * (gs+i2));
    };
};

struct MemSpan {
    bidx_t start[3] = {0, 0, 0};  //!< starting index for the region of interest
    bidx_t end[3]   = {0, 0, 0};  //!< ending index for the region of interest

    explicit MemSpan() = delete;
    explicit MemSpan(const bidx_t in_start, const bidx_t in_end) noexcept;
    explicit MemSpan(const bidx_t in_start[3], const bidx_t in_end[3]) noexcept;
};

/**
 * @brief Translate the ghost limits from one block to another
 * 
 * if the current index is in a ghost region (the front one or the back one), returns the limit of the new ghost region
 * if the current index is in the block, return the index scaled to the new length of the block
 * 
 * @warning if the new index is NOT an integer, we ceil it! (sometimes we use an index of 1 because we want to visit the index 0, we should preserve that)
 * 
 * @param c_id current index in the current reference
 * @param c_core the number of points inside the current block (typically M_N)
 * @param n_front the number of ghost points in front of the block in the new block
 * @param n_core the number of points inside the new block
 * @param n_back the number of ghost points at the back of the block in the new block
 * @return bidx_t the new ghost index
 */
static bidx_t TranslateBlockLimits(const bidx_t c_id, const bidx_t c_core,
                                   const bidx_t n_front, const bidx_t n_core, const bidx_t n_back) {
    // if ratio > 1, we scale up, if ratio < 1, we scale down
    const real_t ratio = static_cast<real_t>(c_core) / static_cast<real_t>(n_core);
    // get where we are (in front, in the block or at the back)
    const bidx_t b = (c_id + c_core);
    const bidx_t c = (b / c_core) + static_cast<bidx_t>(c_id > c_core);
    // compute every possibility
    const bidx_t res[4] = {-n_front, static_cast<bidx_t>(ceil(c_id / ratio)), n_core, n_core + n_back};
    // check that we didn't screw up the indexes
    // m_assert(!(c == 1 && !m_fequal(static_cast<real_t>(c_id) / ratio, res[1])), "we cannot translate the id %d from a core of %d to a core of %d: %f == %d", c_id, c_core, n_core, static_cast<real_t>(c_id) / ratio, res[1]);
    // return the correct choice
    return res[c];
}


#endif