#ifndef SRC_CORE_MEMDATA_HPP_
#define SRC_CORE_MEMDATA_HPP_

#include <type_traits>

#include "core/macros.hpp"
#include "core/memspan.hpp"
#include "core/types.hpp"

// defines a function that takes 3 arguments and returns an offset
// using accessor_t = std::function<bidx_t(const bidx_t, const bidx_t, const bidx_t)>;

/**
 * @brief defines a raw pointer that can be allocated and free'ed
 * 
 */
struct MemPtr {
#ifndef NDEBUG
    bool is_allocated = false;
#endif
    size_t  size = 0;
    real_t* ptr  = nullptr;

    explicit MemPtr() noexcept {};
    ~MemPtr();

    void Allocate(const size_t element_size);
    void Free();

    void MemSetZero();
};

/**
 * @brief defines a __restrict memory region, typically used for the blocks
 * 
 * @warning Do NOT use this class directly, use MemData and ConstMemData instead
 * 
 * @tparam T 
 */
template <typename T>
class RestrictData {
    //--------------------------------------------------------------------------
   private:
    T* __restrict const data_;
    bidx_t const stride_[2]; //!< 1D strides of the data: fastest rotating one and normal one

    //--------------------------------------------------------------------------
    // every class that is like me but with a const type is my friend
    // this allows me to have a move constructor on myself + on a class that is const myself
    friend class RestrictData<const T>;

    //--------------------------------------------------------------------------
   public:
    explicit RestrictData() = delete;
    explicit RestrictData(const MemPtr* ptr, const MemLayout* layout, const lda_t ida = 0) noexcept
        : stride_{layout->stride[0], layout->stride[1]},
          data_(ptr->ptr + layout->offset(ida)) {
        // do alignement asserts
        m_assert(m_isaligned(data_), "the value of data must be aligned!");
        m_assert(ptr->size >= layout->n_elem, "the size of the pointer = %ld must be >= the layout nelem = %ld", ptr->size, layout->n_elem);
    };
    explicit RestrictData(const T* __restrict ptr, const MemLayout* layout, const lda_t ida = 0) noexcept
        : stride_{layout->stride[0], layout->stride[1]},
          data_(ptr + layout->offset(ida)) { m_assert(m_isaligned(data_), "the value of data must be aligned!"); };

    // allow creation of a nullData
    explicit RestrictData(const std::nullptr_t) : data_{nullptr}, stride_{0, 0} {};

    // enable the construction form another MemData<T> (this is a copy construct)
    // note 1: not explicit as it's convenient to not explicity call the constructor
    // note 2: will not compile if T2 is not T or const T (cfr the friendship above)
    template <typename T2>
    RestrictData(const RestrictData<T2>& other) : data_(other.data_), stride_{other.stride_[0], other.stride_[1]} {};

    template <typename T2>
    RestrictData(const RestrictData<T2>& other, const bidx_t i0, const bidx_t i1, const bidx_t i2) noexcept
        : data_(other.ptr(i0, i1, i2)), stride_{other.stride_[0], other.stride_[1]} {};

    // same but with a custom user-given raw address
    // note: will compile whatever T2 is (if it's convertible to T)
    template <typename T2>
    explicit RestrictData(const RestrictData<T2>& other, T2* const data) : data_(data), stride_{other.stride_[0], other.stride_[1]} {};

    //--------------------------------------------------------------------------
    // this is convenient and a bit ugly
    [[nodiscard]] inline bool is_null() const { return data_ == nullptr; };

    //--------------------------------------------------------------------------
    RestrictData<T> shift(const bidx_t i0, const bidx_t i1, const bidx_t i2) const noexcept{
        return RestrictData<T>(ptr(i0,i1,i2));
    }

    //--------------------------------------------------------------------------
    // force the inline! hopefully will do it
    /**
     * @brief returns the offset of a point (i0, i1, i2)
     */
    __attribute__((always_inline)) inline bidx_t offset(const bidx_t i0, const bidx_t i1, const bidx_t i2) const noexcept {
        return (i0 + stride_[0] * (i1 + stride_[1] * i2));
    };

    /**
     * @brief return the data located at (i0,i1,i2) position wrt a reference point given by offset_ref
     * 
     */
    __attribute__((always_inline)) inline T& operator()(const bidx_t i0, const bidx_t i1, const bidx_t i2, const bidx_t offset_ref = 0) const noexcept {
        return data_[offset_ref + i0 + stride_[0] * (i1 + stride_[1] * i2)];
    }

    /**
     * @brief return the data located at (i * ida==0, i * ida==1, i * ida==2) position typically used in stencil operations
     * 
     */
    __attribute__((always_inline)) inline T& Stencil(const bidx_t i, const lda_t ida) const noexcept {
        return data_[(i * (ida == 0)) + stride_[0] * ((i * (ida == 1)) + stride_[1] * (i * (ida == 2)))];
    }

    /**
     * @brief return the pointer to the data located at (i0,i1,i2) position wrt a reference point given by offset_ref
     * 
     */
    __attribute__((always_inline)) inline T* __restrict ptr(const bidx_t i0, const bidx_t i1, const bidx_t i2, const bidx_t offset_ref = 0) const noexcept {
        return data_ + offset_ref + (i0 + stride_[0] * (i1 + stride_[1] * i2));
    }
};

using MemData      = RestrictData<real_t>;
using ConstMemData = RestrictData<const real_t>;

#endif