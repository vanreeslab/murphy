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
using std::string;

/**
 * @brief implements a @ref Block that is used as a leaf for the tree
 * 
 */
class GridBlock : public MemLayout {
   protected:
    bool      lock_;      //!< lock the block, indicating that no refinement/coarsening can happen
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
     * @name Lock management
     * 
     * @{ */
    void lock() { lock_ = true; }
    void unlock() { lock_ = false; }
    bool locked() const { return lock_; }
    /**@} */

    /**
     * @name handle the ghost data pointer
     * @{
     */
    mem_ptr ptr_ghost() const { return ptr_ghost_; }
    void    ptr_ghost(real_p ptr) { ptr_ghost_ = ptr; }
    void    AllocatePtrGhost(const size_t memsize);
    /**@} */

    /**
     * @name datamap access
     * 
     * @{
     */
    // data = memory address of (0,0,0)
    data_ptr data(const Field* fid);
    data_ptr data(const Field* fid, const sid_t ida);
    // data_ptr data(const Field* fid) const;
    // pointer = raw data pointe
    mem_ptr pointer(const Field* fid);
    mem_ptr pointer(const Field* fid, const sid_t ida);
    // mem_ptr pointer(const Field* fid) const;
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


using gridblock_doop_empty_t = void (GridBlock::*)();
using gridblock_doop_field_t = void (GridBlock::*)(Field*);

/**
 * @brief defines a wrapper to call a GridBlock member function as an op_t
 * 
 * @param block the GridBlock
 * @param op the operator of type gop_empty_t
 */
static void CallGridBlockMemFuncEmpty(const qid_t* qid, GridBlock* block, gridblock_doop_empty_t op) {
    (block->*op)();
};

/**
 * @brief defines a wrapper to call a GridBlock member function as an op_t
 * 
 * @param block the GridBlock
 * @param field the considered field
 * @param op the operator of type gop_field_ts
 */
static void CallGridBlockMemFuncField(const qid_t* qid, GridBlock* block, Field* field, gridblock_doop_field_t op) {
    (block->*op)(field);
};


#endif  // SRC_GRIDBLOCK_HPP_
