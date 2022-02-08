// #ifndef SRC_SUBBLOCK_HPP_
// #define SRC_SUBBLOCK_HPP_

// #include "core/macros.hpp"
// #include "core/memlayout.hpp"
// #include "core/types.hpp"

// /**
//  * @brief Translate the ghost limits from one block to another
//  * 
//  * if the current index is in a ghost region (the front one or the back one), returns the limit of the new ghost region
//  * if the current index is in the block, return the index scaled to the new length of the block
//  * 
//  * @warning if the new index is NOT an integer, we ceil it! (sometimes we use an index of 1 because we want to visit the index 0, we should preserve that)
//  * 
//  * @param c_id current index in the current reference
//  * @param c_core the number of points inside the current block (typically M_N)
//  * @param n_front the number of ghost points in front of the block in the new block
//  * @param n_core the number of points inside the new block
//  * @param n_back the number of ghost points at the back of the block in the new block
//  * @return bidx_t the new ghost index
//  */
// static bidx_t TranslateBlockLimits(const bidx_t c_id, const bidx_t c_core,
//                                    const bidx_t n_front, const bidx_t n_core, const bidx_t n_back) {
//     // if ratio > 1, we scale up, if ratio < 1, we scale down
//     const real_t ratio = static_cast<real_t>(c_core) / static_cast<real_t>(n_core);
//     // get where we are (in front, in the block or at the back)
//     const bidx_t b = (c_id + c_core);
//     const bidx_t c = (b / c_core) + static_cast<bidx_t>(c_id > c_core);
//     // compute every possibility
//     const bidx_t res[4] = {-n_front, static_cast<bidx_t>(ceil(c_id / ratio)), n_core, n_core + n_back};
//     // check that we didn't screw up the indexes
//     // m_assert(!(c == 1 && !m_fequal(static_cast<real_t>(c_id) / ratio, res[1])), "we cannot translate the id %d from a core of %d to a core of %d: %f == %d", c_id, c_core, n_core, static_cast<real_t>(c_id) / ratio, res[1]);
//     // return the correct choice
//     return res[c];
// }

// /**
//  * @brief Implementation of a @ref MemLayout as any part of any 3D memory block
//  * 
//  */
// class SubBlock : public MemLayout {
//    protected:
//     bidx_t gs_[2]    = {0, 0};     //!< the ghostsize: in front and at the back of the block
//     bidx_t stride_   = 0;          //!< the stride
//     bidx_t start_[3] = {0, 0, 0};  //!< starting index for the region of interest
//     bidx_t end_[3]   = {0, 0, 0};  //!< ending index for the region of interest

//    public:
//     // constructors
//     explicit SubBlock() = default;

//     explicit SubBlock(const lid_t gs, const lid_t stride);
//     explicit SubBlock(const lid_t gs[2], const lid_t stride);

//     explicit SubBlock(const lid_t gs, const lid_t stride, const lid_t start[3], const lid_t end[3]);
//     explicit SubBlock(const lid_t gs[2], const lid_t stride, const lid_t start[3], const lid_t end[3]);

//     explicit SubBlock(const lid_t gs, const lid_t stride, const lid_t start, const lid_t end);

//     // destructor
//     virtual ~SubBlock() = default;  //!< destructor needed to guarantee the call to the inheritated constructor

//     void Reset(const lid_t gs, const lid_t stride, const lid_t start, const lid_t end);
//     void Reset(const lid_t gs[2], const lid_t stride, const lid_t start[3], const lid_t end[3]);

//     void Resize(const bidx_t new_gs[2], const bidx_t new_core, SubBlock* new_block);

//     /**
//      * @name Memory Layout Implementation
//      * 
//      * @{ */
//     [[nodiscard]] inline bidx_t gs(const lda_t ida) const { return gs_[ida]; }

//     [[nodiscard]] inline bidx_t stride() const override final { return stride_; }
//     [[nodiscard]] inline bidx_t core() const override final { return stride_ - (gs_[0] + gs_[1]); }
//     [[nodiscard]] inline bidx_t gs() const override final { return gs_[0]; }

//     virtual inline bidx_t start(const lda_t ida) const override { return start_[ida]; }
//     virtual inline bidx_t end(const lda_t ida) const override { return end_[ida]; }
//     /** @} */
// };

// #endif  // SRC_SUBBLOCK_HPP_
