#ifndef SRC_CORE_POINTERS_HPP_
#define SRC_CORE_POINTERS_HPP_

#include "core/macros.hpp"
#include "core/pointers.hpp"
#include "core/types.hpp"
#include "core/memlayout.hpp"

/**
 * @brief A class of pointers, might be owned or not
 * 
 * It solves the problem of the owner information about the pointer
 * 
 * An owned pointer has to be free'd while a non owned pointer is a general pointer that cannot be free'd
 * 
 * a few nice refs:
 * - https://www.artima.com/cppsource/rvalue.html
 * - http://thbecker.net/articles/rvalue_references/section_01.html
 * 
 * @tparam T the type to point to
 */
template <typename T>
class m_ptr {
   protected:
    T*   data_;      //!< the underlying information
    bool is_owned_;  //!< indicate if we own the memory and that we should deallocate

   public:
    //-------------------------------------------------------------------------
    // we apply the rule of three: construtor, copy constructor and copy assignment

    /** @brief default constructor */
    explicit m_ptr() noexcept : data_(nullptr), is_owned_(false){};

    /** @brief copy constructor: copy the adress without ownership */
    m_ptr(const m_ptr<T>& ptr) noexcept : data_(ptr()), is_owned_(false){};

    /** @brief copy assignment: copy the address, without ownership */
    m_ptr<T>& operator=(const m_ptr<T>& other) {
        data_     = other();
        is_owned_ = false;
        return *this;
    };

    //-------------------------------------------------------------------------
    // additional usefull constructors

    /** @brief nullptr constructor: any pointer can be initialized to nullptr, no ownership */
    m_ptr(const std::nullptr_t& ptr) noexcept : data_(ptr), is_owned_(false){};

    /** @brief pointer constructor: copy the adress from ptr, no ownership */
    m_ptr(T* ptr) noexcept : data_(ptr), is_owned_(false){};

    /** @brief conversion constructor: copy the adress from ptr, no ownership, fails if not conversion exists */
    template <typename T2>
    m_ptr(const m_ptr<T2>& ptr) noexcept : data_(ptr()), is_owned_(false){};

    //-------------------------------------------------------------------------
    /** @brief destructor, free the pointer if owned */
#ifndef NDEBUG
    ~m_ptr() {
        m_assert(!is_owned_, "attempt to kill a pointer that still owns memory");
    };
#endif

    //-------------------------------------------------------------------------
    /** @brief allocate a new pointer and forward the arguments to the constructor */
    template <typename... U>
    void Alloc(U&&... u) {
        data_     = new T(std::forward<U>(u)...);
        is_owned_ = true;
    };

    /** @brief free the allocated memory if we own it */
    void Free() {
        m_assert(is_owned_, "attempt to free memory but nothing is owned");
        delete (data_);
        data_     = nullptr;
        is_owned_ = false;
    }

    //-------------------------------------------------------------------------
    /** @brief return true if the pointer is nullptr */
    bool IsEmpty() const { return data_ == nullptr; };

    /** @brief return true if it owns the data */
    bool IsOwned() const { return is_owned_; };

    //-------------------------------------------------------------------------
    /** @brief return the underlying pointer */
    T* operator()() const { return data_; }

    /** @brief operator *, return the associated object */
    T& operator*() const { return *data_; }

    /** @brief operator ->, return the contained pointer */
    T* operator->() const { return data_; }
};

/** @brief data pointer type = root of the data, i.e. the adress of (0,0,0)
 * 
 * It relies on a mem_ptr, created first
 */
class data_ptr : public m_ptr<real_t> {
   public:
    using m_ptr<real_t>::m_ptr;       // inheritates the constructors
    using m_ptr<real_t>::operator();  // inheritates the operator()

    real_t* Write(const bidx_t i0 = 0, const bidx_t i1 = 0, const bidx_t i2 = 0, const lda_t ida = 0, const bidx_t stride = M_STRIDE) const;
    real_t* Write(const bidx_t i0, const bidx_t i1, const bidx_t i2, const lda_t ida, const m_ptr<const MemLayout>& layout) const;
    real_t* Write(const lda_t ida, const m_ptr<const MemLayout>& layout) const;
    real_t* Write(const m_ptr<const MemLayout>& layout, const lda_t ida = 0) const;

    const real_t* Read(const bidx_t i0 = 0, const bidx_t i1 = 0, const bidx_t i2 = 0, const lda_t ida = 0, const bidx_t stride = M_STRIDE) const;
    const real_t* Read(const bidx_t i0, const bidx_t i1, const bidx_t i2, const lda_t ida, const m_ptr<const MemLayout>& layout) const;
    const real_t* Read(const lda_t ida, const m_ptr<const MemLayout>& layout) const;
};

/** @brief data pointer type = root of the data, i.e. the adress of (0,0,0)
 * 
 * It relies on a mem_ptr, created first
*/
class const_data_ptr : public m_ptr<const real_t> {
   public:
    using m_ptr<const real_t>::m_ptr;       // inheritates the constructor
    using m_ptr<const real_t>::operator();  // inheritates the operators

    /** @brief build a const_data_ptr from a data_ptr: copy */
    const_data_ptr(const data_ptr& ptr) : m_ptr<const real_t>(ptr){};

    const real_t* Read(const bidx_t i0 = 0, const bidx_t i1 = 0, const bidx_t i2 = 0, const lda_t ida = 0, const bidx_t stride = M_STRIDE) const;
    const real_t* Read(const bidx_t i0, const bidx_t i1, const bidx_t i2, const lda_t ida, const m_ptr<const MemLayout>& layout) const;
    const real_t* Read(const lda_t ida, const m_ptr<const MemLayout>& layout) const;
};

/** @brief pointer to a 3d memory chunck
 * 
 * this is the home of the data, where its allocation (or not is handled)
 * In particular, we make sure that 
 *
 */
class mem_ptr : public m_ptr<real_t> {
    bidx_t stride[2];   //!< array of stride, in the main direction and then in the two other directions
    size_t size   = 0;  //!< store the size of the memory allocated
    size_t offset = 0;  //!< store the offset needed to return a (const_)data_ptr

   public:
    using m_ptr<real_t>::m_ptr;       // inheritates the constructor:
    using m_ptr<real_t>::operator();  // inheritates the operators:

    // define the destructor to free the mem, otherwise it's never done
#ifndef NDEBUG
    ~mem_ptr() {
        m_assert(!is_owned_, "attempt to kill a m_ptr that still own data, call Free() first");
    }
#endif

    /** @brief allocate an owned array of size */
    void Calloc(const size_t size) {
        m_assert(!IsOwned(), "we cannot have some memory already");
        //-------------------------------------------------------------------------
        // allocate the new memory
        data_     = reinterpret_cast<real_t*>(m_calloc(size * sizeof(real_t)));
        is_owned_ = true;
        //-------------------------------------------------------------------------
    }

    /** @brief free the allocated memory if we own it */
    void Free() {
        m_assert(IsOwned(), "attempt to free memory but nothing is owned");
        m_free(data_);
        data_     = nullptr;
        is_owned_ = false;
    }

    data_ptr operator()(const lda_t ida, const bidx_t gs, const bidx_t stride) const;
    data_ptr operator()(const lda_t ida, const m_ptr<const MemLayout>& layout) const;
    mem_ptr  shift_dim(const lda_t ida, const m_ptr<const MemLayout>& layout) const;
};
/**@}*/

// put the assert if we are not in non-debug mode
// #ifdef NDEBUG
constexpr bidx_t m_idx(const bidx_t i0, const bidx_t i1, const bidx_t i2, const bidx_t ida = 0, const bidx_t stride = M_STRIDE) {
    // m_assert((stride*stride*stride) < std::numeric_limits<bidx_t>::max(),"The whole block indexing must be contained in bidx_t");
    bidx_t offset = i0 + stride * (i1 + stride * (i2 + stride * ida));
    return offset;
}
// constexpr bidx_t m_idx(const bidx_t i0, const bidx_t i1, const bidx_t i2, const bidx_t ida, const MemLayout* block) {
//     const bidx_t stride = block->stride();
//     bidx_t       offset = i0 + stride * (i1 + stride * (i2 + stride * ida));
//     return offset;
// }
// #else
// const bidx_t m_idx(const bidx_t i0, const bidx_t i1, const bidx_t i2, const bidx_t ida = 0, const bidx_t stride = M_STRIDE) {
//     m_assert(0 <= stride && stride <= M_STRIDE, "the stride=%d must be between 0 and M_STRIDE", stride);
//     bidx_t offset = i0 + stride * (i1 + stride * (i2 + stride * ida));
//     return offset;
// }
// #endif

// constexpr size_t m_idx(const bidx_t i0, const bidx_t i1, const bidx_t i2, const bidx_t ida, m_ptr<const MemLayout> layout) {
//     const size_t stride = layout->stride();
//     return (size_t)(i0 + stride * (i1 + stride * (i2 + stride * ida)));
// }

#endif