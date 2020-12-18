#include "pointers.hpp"

//=============================================================================
/**
 * @brief return a write access to the data
 * 
 * @param i0 the index in dimension 0 = X
 * @param i1 the index in dimension 1 = Y
 * @param i2 the index in dimension 2 = Z
 * @param ida the dimension
 * @param stride the stride
 * @return real_t* 
 */
real_t* data_ptr::Write(const bidx_t i0, const bidx_t i1, const bidx_t i2, const lda_t ida, const bidx_t stride) const {
    m_assert(0 <= stride && stride <= M_STRIDE, "the stride = %d is wrong", stride);
    //-------------------------------------------------------------------------
    const size_t offset = i0 + stride * (i1 + stride * (i2 + stride * ida));
    return (*this)() + offset;
    //-------------------------------------------------------------------------
}

/**
 * @brief @brief return a write access to the data
 * 
 * see data_ptr::Write(const bidx_t i0, const bidx_t i1, const bidx_t i2, const lda_t ida, const bidx_t stride)
 * 
 * @param i0 the index in dimension 0 = X
 * @param i1 the index in dimension 1 = Y
 * @param i2 the index in dimension 2 = Z
 * @param ida the dimension
 * @param layout the memory layout
 * @return real_t* 
 */
real_t* data_ptr::Write(const bidx_t i0, const bidx_t i1, const bidx_t i2, const lda_t ida, m_ptr<const MemLayout> layout) const {
    //-------------------------------------------------------------------------
    return this->Write(i0, i1, i2, ida, layout()->stride());
    //-------------------------------------------------------------------------
}

/**
 * @brief  return a write access to the data in (0,0,0)
 * 
 * @param ida the given dimension
 * @param layout the memory layout
 * @return const real_t* 
 */
real_t* data_ptr::Write(const lda_t ida, m_ptr<const MemLayout> layout) const {
    m_assert(0 <= layout()->stride() && layout()->stride() <= M_STRIDE, "the stride = %d is wrong", layout()->stride());
    //-------------------------------------------------------------------------
    const size_t stride = layout()->stride();
    const size_t offset = stride * stride * stride * ida;
    return (*this)() + offset;
    //-------------------------------------------------------------------------
}

/**
 * @brief return a read-only access to the data
 * 
 * if i0 = 0, i1 = 0, i2 = 0, and ida = 0, the last argument is discarded
 * 
 * @param i0 the index in dimension 0 = X
 * @param i1 the index in dimension 1 = Y
 * @param i2 the index in dimension 2 = Z
 * @param ida the dimension
 * @param stride the stride
 * @return real_t* 
 */
const real_t* data_ptr::Read(const bidx_t i0, const bidx_t i1, const bidx_t i2, const lda_t ida, const bidx_t stride) const {
    m_assert(0 <= stride && stride <= M_STRIDE, "the stride = %d is wrong", stride);
    //-------------------------------------------------------------------------
    const size_t offset = i0 + stride * (i1 + stride * (i2 + stride * ida));
    return (*this)() + offset;
    //-------------------------------------------------------------------------
}

/**
 * @brief @brief return a read-only access to the data
 * 
 * see data_ptr::Read(const bidx_t i0, const bidx_t i1, const bidx_t i2, const lda_t ida, const bidx_t stride)
 * 
 * @param i0 the index in dimension 0 = X
 * @param i1 the index in dimension 1 = Y
 * @param i2 the index in dimension 2 = Z
 * @param ida the dimension
 * @param layout the memory layout
 * @return real_t* 
 */
const real_t* data_ptr::Read(const bidx_t i0, const bidx_t i1, const bidx_t i2, const lda_t ida, m_ptr<const MemLayout> layout) const {
    //-------------------------------------------------------------------------
    return this->Read(i0, i1, i2, ida, layout()->stride());
    //-------------------------------------------------------------------------
}

/**
 * @brief  return a read-only access to the data in (0,0,0)
 * 
 * @param ida the given dimension
 * @param layout the memory layout
 * @return const real_t* 
 */
const real_t* data_ptr::Read(const lda_t ida, m_ptr<const MemLayout> layout) const {
    m_assert(0 <= layout()->stride() && layout()->stride() <= M_STRIDE, "the stride = %d is wrong", layout()->stride());
    //-------------------------------------------------------------------------
    const size_t stride = layout()->stride();
    const size_t offset = stride * stride * stride * ida;
    return (*this)() + offset;
    //-------------------------------------------------------------------------
}

//=============================================================================
/**
 * @brief return a read-only access to the data
 * 
 * if i0 = 0, i1 = 0, i2 = 0, and ida = 0, the last argument is discarded
 * 
 * @param i0 the index in dimension 0 = X
 * @param i1 the index in dimension 1 = Y
 * @param i2 the index in dimension 2 = Z
 * @param ida the dimension
 * @param layout the memory layout
 * @return real_t* 
 */
const real_t* const_data_ptr::Read(const bidx_t i0, const bidx_t i1, const bidx_t i2, const lda_t ida, const bidx_t stride) const {
    m_assert(0 <= stride && stride <= M_STRIDE, "the stride = %d is wrong", stride);
    //-------------------------------------------------------------------------
    const size_t offset = i0 + stride * (i1 + stride * (i2 + stride * ida));
    return (*this)() + offset;
    //-------------------------------------------------------------------------
}

/**
 * @brief  return a read-only access to the data
 * 
 * @param i0 the index in dimension 0 = X
 * @param i1 the index in dimension 1 = Y
 * @param i2 the index in dimension 2 = Z
 * @param ida the dimension
 * @param layout the memory layout
 * @return const real_t* 
 */
const real_t* const_data_ptr::Read(const bidx_t i0, const bidx_t i1, const bidx_t i2, const lda_t ida, m_ptr<const MemLayout> layout) const {
    m_assert(0 <= layout()->stride() && layout()->stride() <= M_STRIDE, "the stride = %d is wrong", layout()->stride());
    //-------------------------------------------------------------------------
    const size_t stride = layout()->stride();
    const size_t offset = i0 + stride * (i1 + stride * (i2 + stride * ida));
    return (*this)() + offset;
    //-------------------------------------------------------------------------
}

/**
 * @brief  return a read-only access to the data in (0,0,0)
 * 
 * @param ida the given dimension
 * @param layout the memory layout
 * @return const real_t* 
 */
const real_t* const_data_ptr::Read(const lda_t ida, m_ptr<const MemLayout> layout) const {
    m_assert(0 <= layout()->stride() && layout()->stride() <= M_STRIDE, "the stride = %d is wrong", layout()->stride());
    //-------------------------------------------------------------------------
    const size_t stride = layout()->stride();
    const size_t offset = stride * stride * stride * ida;
    return (*this)() + offset;
    //-------------------------------------------------------------------------
}

//=============================================================================
/**
 * @brief return a data_ptr that can be used to access the data, the returned data_ptr does not own the data
 * 
 * @param ida the dimension
 * @param gs the ghost size
 * @param stride the stride
 * @return data_ptr 
 */
data_ptr mem_ptr::operator()(const lda_t ida, const bidx_t gs, const bidx_t stride) const {
    m_assert(0 <= gs && gs <= M_GS, "the gs = %d is wrong", gs);
    m_assert(0 <= stride && stride <= M_STRIDE, "the stride = %d is wrong", stride);
    //-------------------------------------------------------------------------
    // get the offset and return a data_ptr to it
    const size_t offset = gs + stride * (gs + stride * (gs + stride * ida));
    return data_ptr(this->m_ptr::operator()() + offset);
    //-------------------------------------------------------------------------
}

/**
 * @brief return a data_ptr that can be used to access the data, the returned data_ptr does not own the data
 * 
 * @param ida the dimension
 * @param layout the layout used to retrieve the dimension (if not ida = 0)
 * @return data_ptr 
 */
data_ptr mem_ptr::operator()(const lda_t ida, m_ptr<const MemLayout> layout) const {
    m_assert(0 <= layout()->gs() && layout()->gs() <= M_GS, "the gs = %d is wrong", layout()->gs());
    m_assert(0 <= layout()->stride() && layout()->stride() <= M_STRIDE, "the stride = %d is wrong", layout()->stride());
    //-------------------------------------------------------------------------
    // get the offset and return a data_ptr to it
    const size_t gs     = layout()->gs();
    const size_t stride = layout()->stride();
    const size_t offset = gs + stride * (gs + stride * (gs + stride * ida));

    real_t* my_ptr = this->m_ptr::operator()() + offset;
    return data_ptr(my_ptr);
    //-------------------------------------------------------------------------
}

/**
 * @brief return the mem_ptr for the given dimension
 * 
 * @param ida the dimension
 * @param layout the layout, only the stride is used here
 * @return mem_ptr 
 */
mem_ptr mem_ptr::shift_dim(const lda_t ida, m_ptr<const MemLayout> layout) const {
    m_assert(0 <= layout()->stride() && layout()->stride() <= M_STRIDE, "the stride = %d is wrong", layout()->stride());
    //-------------------------------------------------------------------------
    // get the offset and return a data_ptr to it
    const size_t stride = layout()->stride();
    const size_t offset = stride * stride * stride * ida;
    return mem_ptr(this->m_ptr::operator()() + offset);
    //-------------------------------------------------------------------------
}