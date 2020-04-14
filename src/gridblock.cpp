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
    // need to delete the fields not deleted yet
    for (datamap_t::iterator iter = data_map_.begin(); iter != data_map_.end(); iter++) {
        m_log("deleting field %s from the block", iter->first.c_str());
        m_free(iter->second);
    }
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
    m_assert(m_isaligned(data), "M_GS and M_N have to be chosen so that (0,0,0) is aligned in memory");
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

void GridBlock::AddField(const qid_t* qid, Field* fid, nullptr_t ctx) {
    m_begin;
    m_assert(ctx == nullptr, "no context is need in this function");
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
    m_end;
}

void GridBlock::DeleteField(const qid_t* qid, Field* fid, nullptr_t ctx) {
    m_begin;
    m_assert(ctx == nullptr, "no context is need in this function");
    //-------------------------------------------------------------------------
    string name = fid->name();
    m_free(data_map_[name]);
    data_map_.erase(name);
    //-------------------------------------------------------------------------
    m_end;
}