#include "gridblock.hpp"

using std::string;

/**
 * @brief constructs a new Block given a 3D length and a position
 * 
 * @param length the length of the current block (unique number)
 * @param xyz the position of the left,bottom corner (x,y,z)
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
    return data_map_[fid->name()] + m_zeroidx(0, this);
#endif
}
/**
 * @brief return (0,0,0) memory position of a field in a given dimension
 * 
 * @param fid the field
 * @param ida the required dimension
 * @return real_p the memory adress of (0,0,0), hence m_idx(i0,i1,i2) MUST be used
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
    return data_map_[fid->name()] + m_zeroidx(ida, this);
#endif
}

const real_p GridBlock::data(const Field* fid) const {
#ifndef NDEBUG
    datamap_t::const_iterator it = data_map_.find(fid->name());
    m_assert(it != data_map_.end(), "the field %s does not exist in this block", fid->name().c_str());
    real_p data = it->second + m_zeroidx(0, this);
    m_assert(m_isaligned(data), "M_GS and M_N have to be chosen so that (0,0,0) is aligned in memory");
    return data;
#else
    return data_map_[fid->name()] + m_zeroidx(0, this);
#endif
}

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

void GridBlock::AddFields(map<string, Field*>* fields) {
    //-------------------------------------------------------------------------
    // remember if I need to free the memory:
    for (auto iter = fields->begin(); iter != fields->end(); iter++) {
        Field* fid  = iter->second;
        AddField(fid);
    }
    //-------------------------------------------------------------------------
}

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

void GridBlock::DeleteFields() {
    // need to delete the fields not deleted yet
    for (datamap_t::iterator iter = data_map_.begin(); iter != data_map_.end(); iter++) {
        m_free(iter->second);
    }
}