#ifndef SRC_CORE_DATA_HPP
#define SRC_CORE_DATA_HPP

#include <functional>

#include "core/macros.hpp"
#include "core/types.hpp"

/**
 * @brief defines a raw pointer that can be allocated and free'ed
 * 
 */
struct MemPtr {
#ifndef NDEBUG
    size_t size         = 0;
    bool   is_allocated = false;
#endif
    real_t* ptr = nullptr;

    explicit MemPtr() noexcept {};
    ~MemPtr();

    void Allocate(const size_t element_size);
    void Free();
};

struct MemLayout {
    size_t stride[2];  //!< the strides of the data: fastest rotating stride in [0], user stride in [1]
    size_t gs[2];      //!< the front ghost sizes: fastest rotating ghost size in [0], user size in [1]
    size_t n_elem;     //!< the total number of elements in the memory layout

    explicit MemLayout() = delete;
    explicit MemLayout(const lda_t n_dim, const bidx_t n_gs_front, const bidx_t n_block, const bidx_t n_gs_back = -1) noexcept;
};

struct MemSpan {
    bidx_t start[3] = {0, 0, 0};  //!< starting index for the region of interest
    bidx_t end[3]   = {0, 0, 0};  //!< ending index for the region of interest

    explicit MemSpan() = delete;
    explicit MemSpan(const bidx_t in_start, const bidx_t in_end) noexcept;
    explicit MemSpan(const bidx_t in_start[3], const bidx_t in_end[3]) noexcept;
};

// defines a function that takes 3 arguments and returns an offset
using accessor_t = std::function<bidx_t(const bidx_t, const bidx_t, const bidx_t, const bidx_t)>;

class MemData {
   private:
    real_t* const data_;
    size_t const  stride_[2];

   public:
    explicit MemData(const MemLayout& layout, const MemPtr& ptr) noexcept
        : stride_{layout.stride[0], layout.stride[1]},
          data_(ptr.ptr +
                (layout.gs[0] + layout.stride[0] *
                                    (layout.gs[1] + layout.stride[1] * layout.gs[1]))) {
        m_assert(m_isaligned(data_), "the value of data must be aligned!");
    };

    // get now the fancy read, write, etc
    __restrict real_t* RWrite(const lda_t ida = 0) const {
        return data_;
    }
    __restrict const real_t* RRead(const lda_t ida = 0) const {
        return data_;
    }

    real_t* Write(const lda_t ida = 0) const {
        return data_;
    }

    const real_t* Read(const lda_t ida = 0) const {
        return data_;
    }

    [[nodiscard]] __attribute__((always_inline)) inline bidx_t Accessor(const bidx_t i0, const bidx_t i1, const bidx_t i2, const bidx_t ida = 0) const noexcept {
        return i0 + stride_[0] * (i1 + stride_[1] * (i2 + stride_[1] * ida));
    }
};

#endif