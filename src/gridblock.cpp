#include "gridblock.hpp"

using std::string;

/**
 * @brief constructs a new Block given a 3D length and a position
 * 
 * @param length the length of the current block
 * @param xyz the position of the origin, i.e. the left,bottom corner, (x,y,z)
 * @param level the level of the block
 */
GridBlock::GridBlock(const real_t length, const real_t xyz[3], const sid_t level) {
    m_begin;
    //-------------------------------------------------------------------------
    level_ = level;
    for (int id = 0; id < 3; id++) {
        xyz_[id]   = xyz[id];
        hgrid_[id] = length / M_N;
    }
    //-------------------------------------------------------------------------
    m_end;
}

GridBlock::~GridBlock() {
    DeleteFields();
}

/**
 * @brief returns the (aligned!) pointer for write access that corresponds to the first point in the block, i.e. (0,0,0), for the first dimension.
 * You must use either @ref m_sidx, @ref m_midx or @ref m_idx to access any point in the memory
 * 
 * @warning this is not the same pointer as the memory pointers, because the ghost blocks are considered as negative numbers, see @ref MemLayout
 * 
 * @param fid 
 * @return real_p 
 */
real_p GridBlock::data(Field* fid) {
#ifndef NDEBUG
    // check the field validity
    datamap_t::iterator it = data_map_.find(fid->name());
    m_assert(it != data_map_.end(), "the field %s does not exist in this block", fid->name().c_str());
    // check the alignment in memory
    real_p data = it->second + m_zeroidx(0, this);
    m_assert(m_isaligned(data), "M_GS and M_N have to be chosen so that (0,0,0) is aligned in memory");
    return data;
#else
    return data_map_.at(fid->name()) + m_zeroidx(0, this);
#endif
}
/**
 * @brief returns the (aligned!) pointer for write access that corresponds to the first point in the block, i.e. (0,0,0), for the given dimension.
 * You must use either @ref m_sidx, @ref m_midx or @ref m_idx to access any point in the memory
 * 
 * @warning this is not the same pointer as the memory pointers, because the ghost blocks are considered as negative numbers, see @ref MemLayout
 * 
 * @param fid the field
 * @param ida the required dimension
 * @return real_p the memory adress, we ensure its alignement
 */
real_p GridBlock::data(const Field* fid, const sid_t ida) {
#ifndef NDEBUG
    // check the field validity
    datamap_t::iterator it = data_map_.find(fid->name());
    m_assert(it != data_map_.end(), "the field %s does not exist in this block", fid->name().c_str());
    // check the alignment in memory
    real_p data = it->second + m_zeroidx(ida, this);
    m_assert(m_isaligned(data), "M_GS = %d and M_N = %d have to be chosen so that (0,0,0) is aligned in memory: ida = %d -> o", M_GS, M_N, ida);
    return data;
#else
    return data_map_.at(fid->name()) + m_zeroidx(ida, this);
#endif
}

/**
 * @brief  returns the constant (aligned!) pointer for read access that corresponds to the first point in the block, i.e. (0,0,0), for the given dimension.
 * You must use either @ref m_sidx, @ref m_midx or @ref m_idx to access any point in the memory
 * 
 * @warning this is not the same pointer as the memory pointers, because the ghost blocks are considered as negative numbers, see @ref MemLayout
 * 
 * @param fid 
 * @return const real_p the pointer is const
 */
real_p GridBlock::data(const Field* fid) const {
#ifndef NDEBUG
    datamap_t::const_iterator it = data_map_.find(fid->name());
    m_assert(it != data_map_.end(), "the field %s does not exist in this block", fid->name().c_str());
    real_p data = it->second + m_zeroidx(0, this);
    m_assert(m_isaligned(data), "M_GS and M_N have to be chosen so that (0,0,0) is aligned in memory");
    return data;
#else
    return data_map_.at(fid->name()) + m_zeroidx(0, this);
#endif
}

/**
 * @brief adds a field to the block if it doesn't exist already
 * 
 * @param fid 
 */
void GridBlock::AddField(Field* fid) {
    //-------------------------------------------------------------------------
    string name = fid->name();
    // try to find the field
    datamap_t::iterator it = data_map_.find(name);
    // if not found, create it
    if (it == data_map_.end()) {
        m_verb("adding field %s to the block", name.c_str());
        data_map_[name] = (real_p)m_calloc(m_blockmemsize(fid->lda()) * sizeof(real_t));
    } else {
        m_verb("field %s already in the block", name.c_str());
    }
    //-------------------------------------------------------------------------
}

/**
 * @brief add all the fields contained in the map to the current block, if they do not exist already
 * 
 * @param fields 
 */
void GridBlock::AddFields(map<string, Field*>* fields) {
    //-------------------------------------------------------------------------
    // remember if I need to free the memory:
    for (auto iter = fields->begin(); iter != fields->end(); iter++) {
        Field* fid  = iter->second;
        AddField(fid);
    }
    //-------------------------------------------------------------------------
}

/**
 * @brief remove the field from the current block if it exists
 * 
 * @param fid 
 */
void GridBlock::DeleteField(Field* fid) {
    //-------------------------------------------------------------------------
    string name = fid->name();
    // try to find the field
    datamap_t::iterator it = data_map_.find(name);
    // if not found, create it
    if (it != data_map_.end()) {
        m_verb("deleting field %s to the block", name.c_str());
        m_free(data_map_[name]);
        data_map_.erase(name);
    } else {
        m_verb("no field %s in the block", name.c_str());
    }
    //-------------------------------------------------------------------------
}

/**
 * @brief deallocate all the remaining fields in the current block
 */
void GridBlock::DeleteFields() {
    // need to delete the fields not deleted yet
    for (datamap_t::iterator iter = data_map_.begin(); iter != data_map_.end(); iter++) {
        m_free(iter->second);
    }
}
