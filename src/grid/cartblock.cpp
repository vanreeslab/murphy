#include "cartblock.hpp"

using std::string;

/**
 * @brief constructs a new Cartesian block given a 3D length and a position
 * 
 * @param length the length of the current block
 * @param xyz the position of the origin, i.e. the left,bottom corner, (x,y,z)
 * @param level the level of the block
 */
CartBlock::CartBlock(const real_t length, const real_t xyz[3], const level_t level) {
    m_begin;
    //-------------------------------------------------------------------------
    level_ = level;
    for (lda_t id = 0; id < 3; id++) {
        xyz_[id]   = xyz[id];
        hgrid_[id] = CartBlockHGrid(length);// length / (M_N - 1);
    }
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief returns the data_ptr corresponding to a field
 * 
 * @param fid the field
 * @param ida the required dimension
 * @return data_ptr the data_ptr corresponding to the point in the block, i.e. (0,0,0), for the given dimension.
 */
data_ptr CartBlock::data(const m_ptr<const Field>& fid, const lda_t ida) const {
    //-------------------------------------------------------------------------
#ifndef NDEBUG
    // check the field validity
    auto it = mem_map_.find(fid->name());
    m_assert(it != mem_map_.end(), "the field \"%s\" does not exist in this block", fid->name().c_str());
    data_ptr data = it->second(ida, this);
    return data;
#else
    return mem_map_.at(fid->name())(ida, this);
#endif
    //-------------------------------------------------------------------------
}

// /**
//  * @brief  returns the (aligned!) pointer for write access that corresponds to the first point in the block, i.e. (0,0,0), for the first dimension.
//  * You must use either @ref m_sidx, @ref m_midx or @ref m_idx to access any point in the memory
//  *
//  * @param name name of the field you want to access
//  * @return data_ptr
//  */
// data_ptr CartBlock::data(const std::string name, const lda_t ida) {
//     //-------------------------------------------------------------------------
// #ifndef NDEBUG
//     // check the field validity
//     auto it = mem_map_.find(name);
//     m_assert(it != mem_map_.end(), "the field \"%s\" does not exist in this block", name.c_str());
//     // check the alignment in memory
//     data_ptr data = it->second(ida, this) ; //+ m_zeroidx(ida, this);
//     // m_assert(m_isaligned(data), "M_GS = %d and M_N = %d have to be chosen so that (0,0,0) is aligned in memory: ida = %d -> o", M_GS, M_N, ida);
//     return data;
// #else
//     return mem_map_.at(name)(ida, this); // + m_zeroidx(ida, this);
// #endif
//     //-------------------------------------------------------------------------
// }

/**
 * @brief returns the m_ptr that corresponds to the actual data pointer
 * 
 * @warning do not confuse with @ref data() function
 * 
 * @param fid the field
 * @param ida the required dimension (default = 0)
 * @return mem_ptr the memory pointer
 */
mem_ptr CartBlock::pointer(const m_ptr<const Field>& fid, const lda_t ida) const {
#ifndef NDEBUG
    // check the field validity
    auto it = mem_map_.find(fid->name());
    m_assert(it != mem_map_.end(), "the field \"%s\" does not exist in this block", fid->name().c_str());
    // m_assert(m_isaligned(data), "M_GS = %d and M_N = %d have to be chosen so that (0,0,0) is aligned in memory: ida = %d -> o", M_GS, M_N, ida);
    mem_ptr ptr = it->second.shift_dim(ida, this);
    return ptr;
#else
    return mem_map_.at(fid->name());
#endif
}

/**
 * @brief adds a field to the block if it doesn't exist already
 * 
 * @param fid the pointer to the field to add
 */
void CartBlock::AddField(const m_ptr<const Field>& fid) {
    //-------------------------------------------------------------------------
    string name = fid->name();
    // try to find the field
    auto it = mem_map_.find(name);
    // if not found, create it
    if (it == mem_map_.end()) {
        m_verb("adding field <%s> to the block (dim = %d)", name.c_str(), fid->lda());
        mem_map_[name] = mem_ptr();
        auto it        = mem_map_.find(name);
        it->second.Calloc(CartBlockMemNum(fid->lda()));
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
void CartBlock::AddFields(const std::map<string, m_ptr<Field> >* fields) {
    //-------------------------------------------------------------------------
    // remember if I need to free the memory:
    for (auto iter = fields->cbegin(); iter != fields->cend(); iter++) {
        AddField(iter->second);
    }
    //-------------------------------------------------------------------------
}

/**
 * @brief Check weither the gridblock ownes the field 
 * 
 * @param name The name of the requested field
 */

bool CartBlock::IsFieldOwned(const std::string name) const {
    return (mem_map_.find(name) != mem_map_.end());
}

/**
 * @brief remove the field from the current block if it exists
 * 
 * @param fid the field to remove
 */
void CartBlock::DeleteField(const m_ptr<const Field>& fid) {
    //-------------------------------------------------------------------------
    string name = fid->name();

    // try to find the field
    auto it = mem_map_.find(name);

    // if found, delete it
    if (it != mem_map_.end()) {
        m_verb("deleting field <%s> to the block", name.c_str());
        mem_map_.erase(name);
    } else {
        m_verb("no field <%s> in the block", name.c_str());
    }
    //-------------------------------------------------------------------------
}