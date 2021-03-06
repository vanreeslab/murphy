#ifndef SRC_GRID_CARTBLOCK_HPP_
#define SRC_GRID_CARTBLOCK_HPP_

#include <list>
#include <map>
#include <string>

#include "core/memspan.hpp"
#include "core/memdata.hpp"
#include "core/types.hpp"
#include "grid/field.hpp"
#include "core/data.hpp"

// /**
//  * @brief returns the memory size (in # of elements) of a cartesian block
//  * 
//  * @param ida 
//  * @return size_t 
//  */
// constexpr size_t CartBlockMemNum(const lda_t ida) {
//     size_t stride = (M_N + 2 * M_GS);
//     return ida * stride * stride * stride;
// }

/**
 * @brief returns the grid spacing for a cartblock given a block length
 */
constexpr real_t CartBlockHGrid(const real_t length) {
    return length / (real_t)(M_N);
}

/**
 * @brief Cartesian block of data, contains its basic information (location, level,...) and a map to the field contained on the grid
 * 
 * @warning no ghosting is available for that block, see GridBlock
 * 
 */
class CartBlock {
   protected:
    real_t  length_   = 0.0;              //!< the size (length) of the block
    level_t level_    = -1;               //!< the level of the block
    real_t  xyz_[3]   = {0.0, 0.0, 0.0};  //!< the origin of the block
    real_t  hgrid_[3] = {0.0, 0.0, 0.0};  //!< the grid spacing of the block

    std::map<std::string, MemPtr> mem_map_;  //<! a map of the pointers to the actual data
    std::map<std::string, lambda_expr_t > expr_map_;  //<! a map of the expression that might have been stored

   public:
    explicit CartBlock(const real_t length, const real_t xyz[3], const level_t level) noexcept;
    virtual ~CartBlock();

    /**
     * @brief returns a MemLayout corresponding to the this block
     */
    M_INLINE MemLayout BlockLayout() const {
        //----------------------------------------------------------------------
        return MemLayout(M_LAYOUT_BLOCK, M_GS, M_N);
        //----------------------------------------------------------------------
    }
    /**
     * @brief returns a MemSpan for the interior of this block
     */
    M_INLINE MemSpan BlockSpan() const {
        //----------------------------------------------------------------------
        return MemSpan(0, M_N);
        //----------------------------------------------------------------------
    }

    /**
     * @brief returns the absolute position of a given points (i0, i1, i2)
     */
    M_INLINE void pos(const bidx_t i0, const bidx_t i1, const bidx_t i2, real_t m_pos[3]) const {
        //----------------------------------------------------------------------
        m_pos[0] = i0 * hgrid_[0] + xyz_[0];
        m_pos[1] = i1 * hgrid_[1] + xyz_[1];
        m_pos[2] = i2 * hgrid_[2] + xyz_[2];
        //----------------------------------------------------------------------
    }

    /**
     * @brief return the offset for the partitioning buffer
     * 
     * @warning the offset is measured in number of real_t
     */
    virtual size_t PartitionDataOffset() const { return 0; }

    /**
     * @brief pack (store) the needed information in the partitioner buffer
     * 
     * The data will be sent over to another rank during partitioning and recoverd by PartitionDataUnPack()
     * 
     * @warning the information MUST be cast to a real_t
     */
    virtual void PartitionDataPack(real_t* buff) const {};

    /**
     * @brief unpack the needed information from the partitioner buffer
     * 
     * The data has been stored by PartitionDataPack() and send over to another rank
     * 
     * @warning the information MUST be uncast from a real_t
     */
    virtual void PartitionDataUnPack(const real_t* buff){};

    /**
     * @name CartBlock utility functions
     * 
     * @{ */
    [[nodiscard]] inline level_t       level() const { return level_; }
    [[nodiscard]] inline real_t        length() const { return length_; }
    [[nodiscard]] inline real_t        xyz(const int id) const { return xyz_[id]; }
    [[nodiscard]] inline real_t        hgrid(const int id) const { return hgrid_[id]; }
    [[nodiscard]] inline const real_t* hgrid() const { return hgrid_; }
    [[nodiscard]] inline const real_t* xyz() const { return xyz_; }
    /** @} */

    /**
     * @name datamap access
     * @{
     */
    MemData                 data(const Field* fid, const lda_t ida) const noexcept;
    ConstMemData            ConstData(const Field* fid, const lda_t ida) const noexcept;
    real_t* __restrict RawPointer(const Field* const fid, const lda_t ida) const noexcept;

    Data<const real_t>* ConstDataPtr(const Field* const fid, const lda_t ida) const noexcept;

    /** @} */

    /**
     * @name field management
     * 
     * @{
     */
    void AddField(const Field* fid);
    void DeleteField(const Field* fid);
    void AddFields(const std::map<std::string, Field*>* fields);
    bool IsFieldOwned(const std::string& name) const;
    void SetExpr(const Field* field, const lambda_expr_t expr);
    /** @} */
};

// /**
//  * @brief copy the values of a given CartBlock to the temp array
//  *
//  * @param layout the layout to be copied, everything outside the layout is set to 0
//  * @param data the initial data
//  * @param temp the temporary data (of size @ref CartBlockMemNum(1))
//  */
// inline void ToTempMemory(const MemSpan* span, const ConstMemData* data, const MemPtr* temp) {
//     m_begin;
//     //-------------------------------------------------------------------------
//     // reset the whole temp memory size to 0
//     std::memset(temp.ptr, 0,  temp.size * sizeof(real_t));

//     // const real_t* sdata = data.Read();
//     // real_t*       tdata = temp.Write();

//     // for (bidx_t i2 = layout->start(2); i2 < layout->end(2); ++i2) {
//     //     for (bidx_t i1 = layout->start(1); i1 < layout->end(1); ++i1) {
//     //         bidx_t idx = m_idx(0, i1, i2, 0, layout->stride());
//     //         std::memcpy(tdata + idx, sdata + idx, sizeof(real_t) * (layout->end(0) - layout->start(0)));
//     //     }
//     // }
//     auto op = [=](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void{

//     };
//     for_loop(op,span);
//     //-------------------------------------------------------------------------
//     m_end;
// }

//-------------------------------------------------------------------------
// defines code-wide useful lambda functions
// using lambda_i3_t      = std::function<real_t(const bidx_t i0, const bidx_t i1, const bidx_t i2)>;
// using lambda_i3block_t = std::function<real_t(const bidx_t i0, const bidx_t i1, const bidx_t i2, const CartBlock* const block)>;

#endif  //  SRC_GRID_CARTBLOCK_HPP_