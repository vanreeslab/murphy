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
    size_t size         = 0;
    bool   is_allocated = false;
#endif
    real_t* ptr = nullptr;

    explicit MemPtr() noexcept {};
    ~MemPtr();

    void Allocate(const size_t element_size);
    void Free();
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
    size_t const stride_[2];

    //--------------------------------------------------------------------------
    // every class that is like me but with a const type is my friend
    // this allows me to have a move constructor on myself + on a class that is const myself
    friend class RestrictData<const T>;

    //--------------------------------------------------------------------------
   public:
    explicit RestrictData() = delete;
    explicit RestrictData(const MemPtr& ptr, const MemLayout& layout) noexcept
        : stride_{layout.stride[0], layout.stride[1]},
          data_(ptr.ptr + layout.offset(0, 0, 0)) {
        // do alignement asserts
        m_assert(m_isaligned(data_), "the value of data must be aligned!");
        m_assert(ptr.size >= layout.n_elem, "the size of the pointer = %ld must be >= the layout nelem = %ld", ptr.size, layout.n_elem);
    };

    // allow creation of a nullData
    explicit RestrictData(const std::nullptr_t) : data_{nullptr}, stride_{0, 0} {};

    // enable the construction form another MemData<T> (this is a copy construct)
    // note 1: not explicit as it's convenient to not explicity call the constructor
    // note 2: will not compile if T2 is not T or const T (cfr the friendship above)
    template <typename T2>
    RestrictData(const RestrictData<T2>& other) : data_(other.data_), stride_{other.stride_[0], other.stride_[1]} {};

    // same but with a custom user-given raw address
    // note: will compile whatever T2 is (if it's convertible to T)
    template <typename T2>
    explicit RestrictData(const RestrictData<T2>& other, T2* const data) : data_(data), stride_{other.stride_[0], other.stride_[1]} {};

    //--------------------------------------------------------------------------
    // this is convenient and a bit ugly
    [[nodiscard]] inline bool is_null() { return data_ == nullptr; };

    //--------------------------------------------------------------------------
    // force the inline! hopefully will do it
    __attribute__((always_inline)) inline T& operator()(const bidx_t i0, const bidx_t i1, const bidx_t i2) const noexcept {
        return data_[i0 + stride_[0] * (i1 + stride_[1] * i2)];
    }
    __attribute__((always_inline)) inline T* __restrict ptr(const bidx_t i0, const bidx_t i1, const bidx_t i2) const noexcept {
        return data_ + (i0 + stride_[0] * (i1 + stride_[1] * i2));
    }
};

using MemData      = RestrictData<real_t>;
using ConstMemData = RestrictData<const real_t>;

#endif