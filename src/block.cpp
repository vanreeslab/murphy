#include "block.hpp"

using std::string;

/**
 * @brief constructs a new Block given a 3D length and a position
 * 
 * @param length the length of the current block (x,y,z)
 * @param xyz the position of the left,bottom corner (x,y,z)
 * @param level the level of the block
 */
Block::Block(const real_t length[3], const real_t xyz[3],const sid_t level) {
     m_begin;
    //-------------------------------------------------------------------------
    level_ = level;
    for (int id = 0; id < 3; id++) {
        xyz_[id]   = xyz[id];
        hgrid_[id] = length[id] / M_N;
    }
    //-------------------------------------------------------------------------
    m_end;
};

Block::~Block() {
    // need to delete the fields not deleted yet
    for (datamap_t::iterator iter = data_map_.begin(); iter != data_map_.end(); iter++) {
        m_free(iter->second);
    }
}

void Block::AddField(const qid_t* qid, Field* fid, nullptr_t ctx) {
    m_begin;
    m_assert(ctx == nullptr, "no context is need in this function");
    //-------------------------------------------------------------------------
    string name = fid->name();
    // try to find the field
    datamap_t::iterator it = data_map_.find(name);
    // if not found, create it
    if (it == data_map_.end()) {
        m_verb("adding field %s to the block", name.c_str());
        data_map_[name] = (real_p) m_calloc(m_blockmemsize(fid->lda()) * sizeof(real_t));
    } else {
        m_verb("field %s already in the block", name.c_str());
    }
    //-------------------------------------------------------------------------
    m_end;
}

void Block::DeleteField(const qid_t* qid, Field* fid, nullptr_t ctx) {
    m_begin;
    m_assert(ctx == nullptr, "no context is need in this function");
    //-------------------------------------------------------------------------
    string name = fid->name();
    m_free(data_map_[name]);
    data_map_.erase(name);
    //-------------------------------------------------------------------------
    m_end;
}