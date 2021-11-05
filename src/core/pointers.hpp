// #ifndef SRC_CORE_POINTERS_HPP_
// #define SRC_CORE_POINTERS_HPP_

// #include "core/macros.hpp"
// #include "core/types.hpp"

// /**
//  * @brief A class of pointers, might be owned or not
//  * 
//  * It solves the problem of the owner information about the pointer
//  * 
//  * An owned pointer has to be free'd while a non owned pointer is a general pointer that cannot be free'd
//  * 
//  * a few nice refs:
//  * - https://www.artima.com/cppsource/rvalue.html
//  * - http://thbecker.net/articles/rvalue_references/section_01.html
//  * 
//  * @tparam T the type to point to
//  */
// template <typename T>
// class m_ptr {
//    protected:
//     T*   data_;      //!< the underlying information

//    public:
//     //-------------------------------------------------------------------------
//     // we apply the rule of three: construtor, copy constructor and copy assignment

//     /** @brief default constructor */
//     explicit m_ptr() noexcept : data_(nullptr){};

//     /** @brief copy constructor: copy the adress */
//     m_ptr(const m_ptr<T>& ptr) noexcept : data_(ptr()){};

//     /** @brief copy assignment: copy the address */
//     m_ptr<T>& operator=(const m_ptr<T>& other) {
//         data_     = other();
//         return *this;
//     };

//     //-------------------------------------------------------------------------
//     // additional useful constructors

//     /** @brief nullptr constructor: any pointer can be initialized to nullptr */
//     m_ptr(const std::nullptr_t* ptr) noexcept : data_(ptr){};

//     /** @brief pointer constructor: copy the adress from ptr*/
//     m_ptr(T* ptr) noexcept : data_(ptr){};

//     // /** @brief conversion constructor: copy the adress from ptr fails if not conversion exists */
//     // template <typename T2>
//     // m_ptr(const m_ptr<T2>& ptr) noexcept : data_(ptr()){};


//     //-------------------------------------------------------------------------
//     // /** @brief allocate a new pointer and forward the arguments to the constructor */
//     // template <typename... U>
//     // void Alloc(U*&... u) {
//     //     data_     = new T(std::forward<U>(u)...);
//     // };

//     // /** @brief free the allocated memory if we own it */
//     // void Free() {
//     //     delete (data_);
//     //     data_     = nullptr;
//     // }

//     //-------------------------------------------------------------------------
//     /** @brief return true if the pointer is nullptr */
//     bool IsEmpty() const { return data_ == nullptr; };

//     //-------------------------------------------------------------------------
//     /** @brief return the underlying pointer */
//     T* operator()() const { return data_; }

//     // /** @brief operator *, return the associated object */
//     // T* operator*() const { return *data_; }

//     // /** @brief operator ->, return the contained pointer */
//     // T* operator->() const { return data_; }

//     // /** @brief operator *, return the associated object */
//     // T* operator[](const int idx) const { return data_[idx]; }
// };


// /** @brief data pointer type = root of the data, i.e. the adress of (0,0,0)
//  * 
//  * It relies on a mem_ptr, created first
//  */
// class data_ptr : public m_ptr<real_t> {
//    private: 
//     const bidx_t stride_   = 0;                 // the stride of the data 
//     const bidx_t gs_       = 0;                 // the size of the ghost region 

//    public:
//     explicit data_ptr() = default;
//     explicit data_ptr(const bidx_t stride, const bidx_t gs) : m_ptr<real_t>(), stride_(stride), gs_(gs) {};
//     explicit data_ptr(m_ptr<real_t> ptr, const bidx_t stride, const bidx_t gs) : m_ptr<real_t>(ptr), stride_(stride), gs_(gs) {};
//     explicit data_ptr(const data_ptr* other, const bidx_t stride) : m_ptr<real_t>(other()), stride_(stride), gs_(0) {};
//     // explicit data_ptr(const data_ptr* other, const bidx_t stride, const bidx_t gs) : m_ptr<real_t>(other()), stride_(stride), gs_(gs) {};
    
//     /** @brief copy constructor: copy the adress */
//     data_ptr(const data_ptr* ptr) noexcept : m_ptr<real_t>(ptr()), stride_(ptr.stride()), gs_(ptr.gs()) {};
    
//     /** @brief nullptr constructor: any pointer can be initialized to nullptr */
//     explicit data_ptr(const std::nullptr_t* ptr) noexcept : m_ptr<real_t>(ptr), stride_(0), gs_(0){};
    
//     using m_ptr<real_t>::operator();  // inheritates the operator()

//     real_t* __restrict__ Write(const bidx_t i0 = 0, const bidx_t i1 = 0, const bidx_t i2 = 0, const lda_t ida = 0, const bidx_t stride = M_STRIDE) const noexcept;
//     real_t* __restrict__ Write(const bidx_t i0, const bidx_t i1, const bidx_t i2, const lda_t ida, const m_ptr<const MemLayout>& layout) const noexcept;
//     real_t* __restrict__ Write(const lda_t ida, const m_ptr<const MemLayout>& layout) const noexcept;
//     // real_t* Write(const m_ptr<const MemLayout>& layout, const lda_t ida = 0) const;

//     const real_t* __restrict__ Read(const bidx_t i0 = 0, const bidx_t i1 = 0, const bidx_t i2 = 0, const lda_t ida = 0, const bidx_t stride = M_STRIDE) const noexcept;
//     const real_t* __restrict__ Read(const bidx_t i0, const bidx_t i1, const bidx_t i2, const lda_t ida, const m_ptr<const MemLayout>& layout) const noexcept;
//     const real_t* __restrict__ Read(const lda_t ida, const m_ptr<const MemLayout>& layout) const noexcept;
// };

// /** @brief data pointer type = root of the data, i.e. the adress of (0,0,0)
//  * 
//  * It relies on a mem_ptr, created first
// */
// class const_data_ptr : public m_ptr<const real_t> {
//    private: 
//     const bidx_t stride_;             // the stride of the data 
//     const bidx_t gs_;                 // the size of the ghost region 
   
//    public:
//     explicit const_data_ptr() : m_ptr<const real_t>(), stride_(0), gs_(0) {};
//     explicit const_data_ptr(const bidx_t stride, const bidx_t gs) : m_ptr<const real_t>(), stride_(stride), gs_(gs) {};
//     explicit const_data_ptr(m_ptr<const real_t> ptr, const bidx_t stride, const bidx_t gs) : m_ptr<const real_t>(ptr), stride_(stride), gs_(gs) {};
//     explicit const_data_ptr(const const_data_ptr* other, const bidx_t stride) : m_ptr<const real_t>(other()), stride_(stride), gs_(0) {};
    
//     // Conversion from data_ptr to const_data_ptr
//     const_data_ptr(const data_ptr* ptr) : m_ptr<const real_t>(ptr()), stride_(ptr.stride()), gs_(ptr.gs()){};
    
//     /** @brief nullptr constructor: any pointer can be initialized to nullptr */
//     const_data_ptr(const std::nullptr_t* ptr) noexcept : m_ptr<const real_t>(ptr), stride_(0), gs_(0){};

//     using m_ptr<const real_t>::operator();  // inheritates the operators
    
//     inline bidx_t stride() const {return stride_;};
//     inline bidx_t gs() const {return gs_;};

//     const real_t* __restrict__ Read(const bidx_t i0 = 0, const bidx_t i1 = 0, const bidx_t i2 = 0, const lda_t ida = 0, const bidx_t stride = M_STRIDE) const noexcept;
//     const real_t* __restrict__ Read(const bidx_t i0, const bidx_t i1, const bidx_t i2, const lda_t ida, const m_ptr<const MemLayout>& layout) const noexcept;
//     const real_t* __restrict__ Read(const lda_t ida, const m_ptr<const MemLayout>& layout) const noexcept;
// };


// /** @brief pointer to a 3d memory chunck
//  * 
//  * this is the home of the data, where its allocation (or not is handled)
//  * In particular, we make sure that 
//  *
//  */
// class mem_ptr : public m_ptr<real_t> {
//    public:
//     explicit mem_ptr() : m_ptr<real_t>(){};                                           // inheritates the constructor:
//     explicit mem_ptr(real_t* ptr) : m_ptr<real_t>(ptr){};       // inheritates the constructor:
//    /** @brief nullptr constructor: any pointer can be initialized to nullptr */
//     mem_ptr(const std::nullptr_t* ptr) noexcept : m_ptr<real_t>(ptr){};

//     using m_ptr<real_t>::operator();  // inheritates the operators:

//     // define the destructor to free the mem, otherwise it's never done
// #ifndef NDEBUG
//     ~mem_ptr() {}
// #endif

//     /** @brief allocate an owned array of size */
//     void Calloc(const size_t size) {
//         //-------------------------------------------------------------------------
//         // allocate the new memory
//         data_     = reinterpret_cast<real_t*>(m_calloc(size * sizeof(real_t)));
//         // explicitly touch the memory and set 0.0
//         std::memset(data_, 0, size * sizeof(real_t));
//         //-------------------------------------------------------------------------
//     }

//     /** @brief free the allocated memory if we own it */
//     void Free() {
//         m_free(data_);
//         data_     = nullptr;
//     }

//     data_ptr operator()(const lda_t ida, const bidx_t gs = M_GS, const bidx_t stride = M_STRIDE) const noexcept;
//     data_ptr operator()(const lda_t ida, const m_ptr<const MemLayout>& layout) const noexcept;
//     mem_ptr  shift_dim(const lda_t ida, const m_ptr<const MemLayout>& layout) const noexcept;
// };
// /**@}*/

// // put the assert if we are not in non-debug mode
// // #ifdef NDEBUG

// /**
//  * @brief return the shift in memory while using data_ptr
//  * 
//  * @param i0 the index in the X direction (fastest rotating index)
//  * @param i1 the index in the Y direction
//  * @param i2 the index in the Z direction
//  * @param ida the index of the dimension (default: 0)
//  * @param stride the stride (default: M_STRIDE)
//  * @return constexpr bidx_t 
//  */
// constexpr bidx_t m_idx(const bidx_t i0, const bidx_t i1, const bidx_t i2, const bidx_t ida = 0, const bidx_t stride = M_STRIDE) {
//     bidx_t offset = i0 + stride * (i1 + stride * (i2 + stride * ida));
//     return offset;
// }

// /**
//  * @brief returns an increment of index in one direction
//  * 
//  * This function has been designed to easily handle the stencil computations (for cross shapes)
//  * 
//  * @warning it is not possible to change the dimension with this function
//  * 
//  * @param id the increment of index 
//  * @param ida_de the direction in which the increment happens
//  * @param stride the custom stride (default: M_STRIDE)
//  * @return constexpr bidx_t 
//  */
// constexpr bidx_t m_idx_delta(const bidx_t delta, const lda_t ida_delta, const bidx_t stride = M_STRIDE) {
//     bidx_t i0     = delta * (ida_delta == 0);
//     bidx_t i1     = delta * (ida_delta == 1);
//     bidx_t i2     = delta * (ida_delta == 2);
//     bidx_t offset = i0 + stride * (i1 + stride * i2);
//     return offset;
// }
// // constexpr bidx_t m_idx(const bidx_t i0, const bidx_t i1, const bidx_t i2, const bidx_t ida, const MemLayout* block) {
// //     const bidx_t stride = block->stride();
// //     bidx_t       offset = i0 + stride * (i1 + stride * (i2 + stride * ida));
// //     return offset;
// // }
// // #else
// // const bidx_t m_idx(const bidx_t i0, const bidx_t i1, const bidx_t i2, const bidx_t ida = 0, const bidx_t stride = M_STRIDE) {
// //     m_assert(0 <= stride && stride <= M_STRIDE, "the stride=%d must be between 0 and M_STRIDE", stride);
// //     bidx_t offset = i0 + stride * (i1 + stride * (i2 + stride * ida));
// //     return offset;
// // }
// // #endif

// // constexpr size_t m_idx(const bidx_t i0, const bidx_t i1, const bidx_t i2, const bidx_t ida, m_ptr<const MemLayout> layout) {
// //     const size_t stride = layout->stride();
// //     return (size_t)(i0 + stride * (i1 + stride * (i2 + stride * ida)));
// // }

// #endif
