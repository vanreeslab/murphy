#ifndef SRC_CORE_POINTERS_HPP_
#define SRC_CORE_POINTERS_HPP_

#include "core/types.hpp"
#include "core/macros.hpp"
#include "core/memlayout.hpp"

/**
 * @brief A class of pointers, might be owned or not
 * 
 * It solves two problem:
 * - the difference between mem_ptr and data_ptr by using a strong typing
 * - the owner information about the pointer
 * 
 * An owned pointer is an array and will be allocated at the instantiation
 * A non owned pointer is a general pointer that cannot be free'd
 * 
 * @tparam T the type to point to
 */
template <typename T>
class m_ptr {
   protected:
    T*   data_     = nullptr;  //!< the underlying information
    bool is_owned_ = false;    //!< indicate if we own the memory and that we should deallocate

   public:
    //-------------------------------------------------------------------------
    /** @brief default constructor */
    m_ptr() : data_(nullptr), is_owned_(false){};

    /** @brief copy constructor: copy the pointer and don't own the data */
    template <typename T2>
    m_ptr(const m_ptr<T2>& ptr) : data_(ptr()), is_owned_(false){};

    /** @brief move constructor: deleted */
    template <typename T2>
    m_ptr(m_ptr<T2>&& ptr) : data_(ptr()), is_owned_(false){};

    //-------------------------------------------------------------------------
    /** @brief nullptr constructor */
    m_ptr(std::nullptr_t& ptr) : data_(nullptr), is_owned_(false){};

    /** @brief from pointer constructor: copy the adress from ptr but doesn't own it */
    m_ptr(T* ptr) : data_(ptr), is_owned_(false){};

    /** @brief from size constructor: allocate an owned array of size */
    m_ptr(const size_t size) : data_(reinterpret_cast<T*>(m_calloc(size * sizeof(T)))),
                               is_owned_(true){};

    //-------------------------------------------------------------------------
    /** @brief destructor, free the array if owned */
    ~m_ptr() {
        Free();
    };

    //-------------------------------------------------------------------------
    /** @brief return true if the pointer is nullptr */
    bool IsEmpty() const { return data_ == nullptr; };

    /** @brief return true if it owns the data */
    bool IsOwned() const { return is_owned_; };

    /** @brief return the underlying pointer */
    T* operator()() const { return data_; }

    /** @brief allocate an owned array of size */
    void Calloc(const size_t size) {
        m_assert(!(data_ != nullptr && is_owned_ != true), "we cannot overwrite this pointer, we would loose its content");
        //-------------------------------------------------------------------------
        this->Free();
        // allocate the new memory
        data_     = reinterpret_cast<T*>(m_calloc(size * sizeof(T)));
        is_owned_ = true;
        //-------------------------------------------------------------------------
    }

    /** @brief free the allocated memory if we own it */
    void Free() {
        if (is_owned_) {
            m_free(data_);
            data_ = nullptr;
            is_owned_ = false;
        }
    }
};

/** @brief data pointer type = root of the data, i.e. the adress of (0,0,0) */
class data_ptr : public m_ptr<real_t> {
   public:
    using m_ptr<real_t>::m_ptr;       // inheritates the constructor:
    using m_ptr<real_t>::operator();  // inheritates the operators:

    real_t* Write(const bidx_t i0 = 0, const bidx_t i1 = 0, const bidx_t i2 = 0, const lda_t ida = 0, const bidx_t stride = M_STRIDE) const;
    real_t* Write(const bidx_t i0, const bidx_t i1, const bidx_t i2, const lda_t ida, m_ptr<const MemLayout> layout) const;
    real_t* Write(const lda_t ida, m_ptr<const MemLayout> layout) const;

    const real_t* Read(const bidx_t i0 = 0, const bidx_t i1 = 0, const bidx_t i2 = 0, const lda_t ida = 0, const bidx_t stride = M_STRIDE) const;
    const real_t* Read(const bidx_t i0, const bidx_t i1, const bidx_t i2, const lda_t ida, m_ptr<const MemLayout> layout) const;
    const real_t* Read(const lda_t ida, m_ptr<const MemLayout> layout) const;
};

/** @brief data pointer type = root of the data, i.e. the adress of (0,0,0) */
class const_data_ptr : public m_ptr<const real_t> {
   public:
    using m_ptr<const real_t>::m_ptr;       // inheritates the constructor:
    using m_ptr<const real_t>::operator();  // inheritates the operators:

    /** @brief build a const_data_ptr from a data_ptr */
    const_data_ptr(const data_ptr& ptr) : m_ptr<const real_t>(ptr){};

    const real_t* Read(const bidx_t i0 = 0, const bidx_t i1 = 0, const bidx_t i2 = 0, const lda_t ida = 0, const bidx_t stride = M_STRIDE) const;
    const real_t* Read(const bidx_t i0, const bidx_t i1, const bidx_t i2, const lda_t ida, m_ptr<const MemLayout> layout) const;
    const real_t* Read(const lda_t ida, m_ptr<const MemLayout> layout) const;
};

/** @brief memory pointer type = root of the data, i.e. the adress of (0,0,0) */
class mem_ptr : public m_ptr<real_t> {
   public:
    using m_ptr<real_t>::m_ptr;       // inheritates the constructor:
    using m_ptr<real_t>::operator();  // inheritates the operators:

    data_ptr operator()(const lda_t ida, const bidx_t gs = M_GS, const bidx_t stride = M_STRIDE) const;
    data_ptr operator()(const lda_t ida, m_ptr<const MemLayout> layout) const;

    mem_ptr shift_dim(const lda_t ida, m_ptr<const MemLayout> layout) const;
};

// using data_ptr      = m_ptr<real_t>;        //!< data pointer type = root of the data, i.e. the adress of (0,0,0)
// using mem_ptr       = m_ptr<real_t>;        //!< pointer type = root of the memory allocation
/**@}*/

constexpr bidx_t m_idx(const bidx_t i0, const bidx_t i1, const bidx_t i2, const bidx_t ida = 0, const bidx_t stride = M_STRIDE) {
    bidx_t offset = i0 + stride * (i1 + stride * (i2 + stride * ida));
    return offset;
}

// constexpr size_t m_idx(const bidx_t i0, const bidx_t i1, const bidx_t i2, const bidx_t ida, m_ptr<const MemLayout> layout) {
//     const size_t stride = layout()->stride();
//     return (size_t)(i0 + stride * (i1 + stride * (i2 + stride * ida)));
// }

#endif