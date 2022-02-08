#ifndef SRC_CORE_DATA_HPP_
#define SRC_CORE_DATA_HPP_

/**
 * @brief An abstract class that defines the interface to any data of type <T>
 * 
 * It implements the stencil capability through the shift command
 * It also defines the standard interface to access the data
 * 
 * @tparam T the type of data
 */
template <typename T>
class Data {
   public:
    /**
     * @brief return the data located at (i0,i1,i2) position wrt a reference point given by offset_ref
     * 
     * Must be overwritten by every subclass
     * 
     */
    virtual M_INLINE T at(const bidx_t i0, const bidx_t i1, const bidx_t i2) const noexcept {
        m_assert(false, "this function must be overwritten by the classes and cannot be used this way");
        return 0.0;
    };
};

class LocalData {
    const bidx_t                    offset_[3];  //!< offset at which the stencil is evaluated
    const Data<const real_t>* const data_;       //!< original data sources

   public:
    explicit LocalData() = delete;
    explicit LocalData(const Data<const real_t>* data, const bidx_t i0, const bidx_t i1, const bidx_t i2) noexcept
        : offset_{i0, i1, i2}, data_(data) {}

    /**
     * @brief return the data located at (i * ida==0, i * ida==1, i * ida==2) position typically used in stencil operations
     */
    M_INLINE real_t operator()(const bidx_t i, const lda_t ida) const noexcept {
        //----------------------------------------------------------------------
        return data_->at(offset_[0] + i * (ida == 0),
                         offset_[1] + i * (ida == 1),
                         offset_[2] + i * (ida == 2));
        //----------------------------------------------------------------------
    };
    /**
     * @brief return the data located at the given position, encoded from 0 to 8
     */
    M_INLINE real_t operator()(const bidx_t i0, const bidx_t i1, const bidx_t i2) const noexcept {
        //----------------------------------------------------------------------
        return data_->at(offset_[0] + i0,
                         offset_[1] + i1,
                         offset_[2] + i2);
        //----------------------------------------------------------------------
    };
};

// defined for convenience only
// using Data<const real_t> = Data<const real_t>;

#endif  // SRC_CORE_DATA_HPP_