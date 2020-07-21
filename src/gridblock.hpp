#ifndef SRC_GRIDBLOCK_HPP_
#define SRC_GRIDBLOCK_HPP_

#include <limits>
#include <map>
#include <string>

#include "field.hpp"
#include "memlayout.hpp"
#include "murphy.hpp"
#include "p8est.h"

using std::map;
using std::numeric_limits;

/**
 * @brief implements a @ref Block that is used as a leaf for the tree
 * 
 */
class GridBlock : public MemLayout {
   protected:
    level_t   level_;     //!< the level of the block
    real_t    xyz_[3];    //!< the origin of the block
    real_t    hgrid_[3];  //!< the grid spacing of the block
    datamap_t data_map_;  //<! a map of the pointers to the actual data

    real_t* ptr_ghost_ = nullptr;  //!< a pointer of data that I do not handle but which uniquely associated to me for the ghost computation

   public:
    GridBlock(const real_t length, const real_t xyz[3], const sid_t level);
    ~GridBlock();

    /**
     * @name Memory Layout Implementation
     * 
     * @{ */
    inline lid_t gs() const override { return M_GS; }
    inline lid_t stride() const override { return M_STRIDE; }
    inline lid_t start(const int id) const override { return 0; }
    inline lid_t end(const int id) const override { return M_N; }
    /** @} */

    inline sid_t  level() const { return level_; }
    inline real_t xyz(const int id) const { return xyz_[id]; }
    inline real_t hgrid(const int id) const { return hgrid_[id]; }
    const real_t* hgrid() const { return hgrid_; }
    const real_t* xyz() const { return xyz_; }

    /**
     * @name handle the ghost data pointer
     * @{
     */
    void   AllocatePtrGhost(const size_t memsize);
    void   ptr_ghost(real_p ptr) { ptr_ghost_ = ptr; }
    real_p ptr_ghost() const { return ptr_ghost_; }
    /**@} */

    /**
     * @name datamap access
     * 
     * @{
     */
    // data = memory address of (0,0,0)
    real_p data(Field* fid);
    real_p data(const Field* fid, const sid_t ida);
    real_p data(const Field* fid) const;
    // pointer = raw data pointe
    real_p pointer(Field* fid);
    real_p pointer(const Field* fid, const sid_t ida);
    real_p pointer(const Field* fid) const;
    /** @} */

    /**
     * @name field management
     * 
     * @{
     */
    void AddField(Field* fid);
    void DeleteField(Field* fid);
    void AddFields(map<string, Field*>* fields);
    void DeleteFields();
    /** @} */
};

/**
 * @brief pointer to an member function of the class @ref GridBlock
 */
using gbop_t = void (GridBlock::*)(Field* fid);

#endif  // SRC_GRIDBLOCK_HPP_
