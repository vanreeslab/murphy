#ifndef SRC_GRID_CARTBLOCK_HPP_
#define SRC_GRID_CARTBLOCK_HPP_

#include <list>
#include <map>
#include <string>


// #include "core/pointers.hpp"
#include "core/memspan.hpp"
#include "core/memdata.hpp"
#include "core/types.hpp"
#include "grid/field.hpp"

/**
 * @brief returns the memory size (in # of elements) of a cartesian block
 * 
 * @param ida 
 * @return size_t 
 */
constexpr size_t CartBlockMemNum(const lda_t ida) {
    size_t stride = (M_N + 2 * M_GS);
    return ida * stride * stride * stride;
}

constexpr real_t CartBlockHGrid(const real_t length) {
    // return length / (M_N - 1);
    return length / (real_t)(M_N);
}

/**
 * @brief Cartesian block of data, contains its location, level, a map with fields and a grid spacing
 * 
 * 
 * We consider a cartesian block where the number of point is EVEN
 * A cartesian block is schematically represented as:
 * 
 * ```
 *      M_GS                  M_N                   M_GS
 *   <------>  <------------------------------>   <------>
 *  |         |                                  |         |
 *  x----x----o----o----o----o----o----o----o----x----x----|
 *  |         |                                  |         |
 * ```
 * 
 * @warning two neighbors share the same first/last information!
 * @warning no ghosting is available for that block
 * 
 */
class CartBlock{
   protected:
    level_t level_    = -1;               //!< the level of the block
    real_t  xyz_[3]   = {0.0, 0.0, 0.0};  //!< the origin of the block
    real_t  hgrid_[3] = {0.0, 0.0, 0.0};  //!< the grid spacing of the block

    std::map<std::string, MemPtr> mem_map_;  //<! a map of the pointers to the actual data

   public:
    explicit CartBlock(const real_t length, const real_t xyz[3], const level_t level);

    /**
     * @name Memory Layout functions
     * 
     * @{ */
    inline bidx_t gs() const override final { return M_GS; }
    inline bidx_t core() const override final { return M_N; }
    inline bidx_t stride() const override final { return (M_N + 2 * M_GS); }
    inline bidx_t start(const lda_t ida) const override final { return 0; }
    inline bidx_t end(const lda_t ida) const override final { return M_N; }
    /** @} */

    /**
     * @brief returns the position of a given points (i0, i1, i2)
     */
    inline void pos(const bidx_t i0, const bidx_t i1, const bidx_t i2, real_t m_pos[3]) const {
        m_pos[0] = i0 * hgrid_[0] + xyz_[0];
        m_pos[1] = i1 * hgrid_[1] + xyz_[1];
        m_pos[2] = i2 * hgrid_[2] + xyz_[2];
    }

    /**
     * @name CartBlock utility functions
     * 
     * @{ */
    inline level_t level() const { return level_; }
    inline real_t  xyz(const int id) const { return xyz_[id]; }
    inline real_t  hgrid(const int id) const { return hgrid_[id]; }
    const real_t*  hgrid() const { return hgrid_; }
    const real_t*  xyz() const { return xyz_; }
    /** @} */

    /**
     * @name datamap access
     * @{
     */
    MemData data(const Field* const fid, const lda_t ida = 0) const noexcept;
    // MemPtr  pointer(const Field* const fid, const lda_t ida = 0) const noexcept;
    /** @} */

    /**
     * @name field management
     * 
     * @{
     */
    void AddField(const Field* const fid);
    void DeleteField(const Field* const fid);
    void AddFields(const std::map<std::string, Field*>* fields);
    bool IsFieldOwned(const std::string name) const;
    /** @} */
};

/**
 * @brief copy the values of a given CartBlock to the temp array
 * 
 * @param layout the layout to be copied, everything outside the layout is set to 0
 * @param data the initial data
 * @param temp the temporary data (of size @ref CartBlockMemNum(1))
 */
inline void ToTempMemory(MemLayout* const  layout, const ConstMemData& data, const MemData& temp) {
    m_begin;
    //-------------------------------------------------------------------------
    // reset the memory size to 0
    std::memset(temp(), 0, sizeof(real_t) * CartBlockMemNum(1));

    const real_t* sdata = data.Read();
    real_t*       tdata = temp.Write();

    for (bidx_t i2 = layout->start(2); i2 < layout->end(2); ++i2) {
        for (bidx_t i1 = layout->start(1); i1 < layout->end(1); ++i1) {
            bidx_t idx = m_idx(0, i1, i2, 0, layout->stride());
            std::memcpy(tdata + idx, sdata + idx, sizeof(real_t) * (layout->end(0) - layout->start(0)));
        }
    }
    //-------------------------------------------------------------------------
    m_end;
}

//-------------------------------------------------------------------------
// defines code-wide usefull lambda functions
// using lambda_i3_t      = std::function<real_t(const bidx_t i0, const bidx_t i1, const bidx_t i2)>;
// using lambda_i3block_t = std::function<real_t(const bidx_t i0, const bidx_t i1, const bidx_t i2, const CartBlock* const block)>;

#endif  //  SRC_GRID_CARTBLOCK_HPP_