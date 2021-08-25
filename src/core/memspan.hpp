#ifndef SRC_CORE_MEMSPAN_HPP_
#define SRC_CORE_MEMSPAN_HPP_

#include "core/types.hpp"
#include "core/macros.hpp"

typedef enum m_layout_t {
    M_LAYOUT_BLOCK = 3
} m_layout_type_t;

struct MemLayout {
    bidx_t gs;         //!< number of ghost points in front of the block
    bidx_t shift;      //!< the shift to apply to the fastest rotating dimension to have the 0 aligned
    bidx_t block;      //!< the number of points in the block
    bidx_t stride[2];  //!< the strides of the data: fastest rotating stride in [0], user stride in [1]

    size_t n_elem;  //!< number of elements in the memory layout (= minimal allocation size, can be large and always >= 0)

    explicit MemLayout() = delete;
    explicit MemLayout(const m_layout_t layout, const bidx_t n_gs_front, const bidx_t n_block, const bidx_t n_gs_back = -1) noexcept;

    bidx_t TranslateLimit(const MemLayout* new_layout, const bidx_t id) const noexcept;

    /**
     * @brief returns the offset of the position (0,0,0) of a given dimension from the RAW pointer expressed as a # of elements (>=0 number)
     */
    __attribute__((always_inline)) inline size_t offset(const lda_t ida = 0) const noexcept {
        m_assert(0 <= ida, "ida = %d cannot be <0", ida);
        return shift + gs + stride[0] * (gs + stride[1] * (gs + stride[1] * ida));
    };

    /**
     * @brief returns the offset of the position (i0,i1,i2) from the RAW pointer expressed as a # of elements (>=0 number)
     */
    __attribute__((always_inline)) inline size_t offset(const bidx_t i0, const bidx_t i1, const bidx_t i2) const noexcept {
        m_assert(0 <= (gs + i0) && (gs + i0) < stride[1], "we cannot be <0 or > stride from the user");
        m_assert(0 <= (gs + i1) && (gs + i1) < stride[1], "we cannot be <0 or > stride from the user");
        m_assert(0 <= (gs + i2) && (gs + i2) < stride[1], "we cannot be <0 or > stride from the user");
        return shift + (gs + i0) + stride[0] * ((gs + i1) + stride[1] * (gs + i2));
    };
};

struct MemSpan {
    bidx_t start[3] = {0, 0, 0};  //!< starting index for the region of interest
    bidx_t end[3]   = {0, 0, 0};  //!< ending index for the region of interest

    explicit MemSpan(){};
    explicit MemSpan(const bidx_t in_start, const bidx_t in_end) noexcept;
    explicit MemSpan(const bidx_t in_start[3], const bidx_t in_end[3]) noexcept;
    explicit MemSpan(const MemSpan* old_span, const bidx_t shift[3]) noexcept;
};

// void Translate(const MemLayout* old_layout, const MemSpan* old_span,
//                const MemLayout* new_layout, MemSpan* new_span);

#endif