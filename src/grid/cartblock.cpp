#include "cartblock.hpp"

using std::string;

/**
 * @brief constructs a new Cartesian block given a 3D length and a position
 * 
 * @param length the length of the current block
 * @param xyz the position of the origin, i.e. the left,bottom corner, (x,y,z)
 * @param level the level of the block
 */
CartBlock::CartBlock(const real_t length, const real_t xyz[3], const level_t level) noexcept {
    m_begin;
    //-------------------------------------------------------------------------
    level_  = level;
    length_ = length;

#pragma unroll 3
    for (lda_t id = 0; id < 3; id++) {
        xyz_[id]   = xyz[id];
        hgrid_[id] = CartBlockHGrid(length);  // length / (M_N - 1);
    }
    //--------------------------------------------------------------------------
    m_end;
}

CartBlock::~CartBlock() {
    //--------------------------------------------------------------------------
    // Free the blocks in the mapping
    for (auto it = mem_map_.begin(); it != mem_map_.end(); ++it) {
        it->second.Free();
    }
    //--------------------------------------------------------------------------
}

/**
 * @brief returns the MemData corresponding to a field
 * 
 * @param fid the field
 * @param ida the required dimension
 * @return MemData the MemData corresponding to the origin of the block, i.e. (0,0,0), for the given dimension.
 */
MemData CartBlock::data(const Field* fid, const lda_t ida) const noexcept {
    // warning, from cpp reference
    // Non-throwing functions are permitted to call potentially-throwing functions.
    // Whenever an exception is thrown and the search for a handler encounters the outermost block of a non-throwing function, the function std::terminate or std::unexpected (until C++17) is called
    //-------------------------------------------------------------------------
    MemLayout myself = BlockLayout();
#ifndef NDEBUG
    // check the field validity
    auto it = mem_map_.find(fid->name());
    m_assert(it != mem_map_.end(), "the field \"%s\" does not exist in this block", fid->name().c_str());
    MemData data_out(&it->second, &myself, ida);
#else
    MemData data_out(&mem_map_[fid->name()],&myself, ida);
#endif
    return data_out;
    //-------------------------------------------------------------------------
}
/**
 * @brief returns the ConstMemData corresponding to a field 
 * 
 * @param fid 
 * @param ida 
 * @return ConstMemData the ConstMemData corresponding to the origin of the block, i.e. (0,0,0), for the given dimension.
 */
ConstMemData CartBlock::ConstData(const Field* fid, const lda_t ida) const noexcept {
    // warning, from cpp reference
    // Non-throwing functions are permitted to call potentially-throwing functions.
    // Whenever an exception is thrown and the search for a handler encounters the outermost block of a non-throwing function, the function std::terminate or std::unexpected (until C++17) is called
    //-------------------------------------------------------------------------
    MemLayout myself = BlockLayout();
#ifndef NDEBUG
    // check the field validity
    auto it = mem_map_.find(fid->name());
    m_assert(it != mem_map_.end(), "the field \"%s\" does not exist in this block", fid->name().c_str());
    ConstMemData data_out(&it->second, &myself, ida);
#else
    ConstMemData data_out(&mem_map_[fid->name()], &myself, ida);
#endif
    return data_out;
    //-------------------------------------------------------------------------
}


/**
 * @brief returns the raw pointer corresponding to a field
 * 
 * @param fid the field
 * @param ida the required dimension
 * 
 * @return real_t* the raw pointer corresponding to the given field. 
 *         useful for partitionning and for the IO
 *
 * @warning Do not use this function unless you are sure about what you are doing 
 * 
 */
real_t* __restrict CartBlock::RawPointer(const Field* fid, const lda_t ida) const noexcept { 
    // warning, from cpp reference
    // Non-throwing functions are permitted to call potentially-throwing functions.
    // Whenever an exception is thrown and the search for a handler encounters the outermost block of a non-throwing function, the function std::terminate or std::unexpected (until C++17) is called
    //-------------------------------------------------------------------------
#ifndef NDEBUG
    // check the field validity
    auto it = mem_map_.find(fid->name());
    m_assert(it != mem_map_.end(), "the field \"%s\" does not exist in this block", fid->name().c_str());
#endif
    return mem_map_.at(fid->name()).ptr + BlockLayout().n_elem*ida;
}


/**
 * @brief adds a field to the block if it doesn't exist already
 * 
 * @param fid the pointer to the field to add
 */
void CartBlock::AddField(const Field* fid) {
    //-------------------------------------------------------------------------
    string name = fid->name();
    // try to find the field
    auto it = mem_map_.find(name);

    // if not found, create it
    if (it == mem_map_.end()) [[likely]] {
        // create an empty pointer, a copy is happening here!
        mem_map_[name] = MemPtr();

        // allocate the pointer, the one in the map, not the one created
        MemLayout myself = BlockLayout();
        auto      it = mem_map_.find(name);
        it->second.Allocate(myself.n_elem * fid->lda());

        m_verb("adding field <%s> to the block (dim = %d)", name.c_str(), fid->lda());
    } else {
        m_verb("field <%s> already in the block (dim=%d)", name.c_str(), fid->lda());
    }
    //-------------------------------------------------------------------------
}

/**
 * @brief recursivelly call @ref AddField() on everyfield in the mapping
 * 
 * @param fields 
 */
void CartBlock::AddFields(const std::map<string, Field*>* fields) {
    //-------------------------------------------------------------------------
    // remember if I need to free the memory:
    for (const auto fid : (*fields) ) {
        AddField(fid.second);
    }
    //-------------------------------------------------------------------------
}

/**
 * @brief Check weither the gridblock ownes the field 
 * 
 * @param name The name of the requested field
 */

bool CartBlock::IsFieldOwned(const string& name) const {
    return (mem_map_.find(name) != mem_map_.end());
}

/**
 * @brief remove the field from the current block if it exists
 * 
 * @param fid the field to remove
 */
void CartBlock::DeleteField(const Field* fid) {
    //-------------------------------------------------------------------------
    string name = fid->name();

    // try to find the field
    auto it = mem_map_.find(name);

    // if found, delete it
    if (it != mem_map_.end()) {
        it->second.Free();
        m_verb("deleting field <%s> to the block", name.c_str());
        mem_map_.erase(name);
    } else {
        m_verb("no field <%s> in the block", name.c_str());
    }
    //-------------------------------------------------------------------------
}