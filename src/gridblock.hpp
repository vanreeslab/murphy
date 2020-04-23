#ifndef SRC_BLOCK_HPP_
#define SRC_BLOCK_HPP_

#include <limits>
#include <map>

#include "field.hpp"
#include "memlayout.hpp"
#include "murphy.hpp"
#include "p8est.h"

using std::map;
using std::numeric_limits;

class GridBlock : public MemLayout {
   protected:
    sid_t  level_;
    real_t xyz_[3];
    real_t hgrid_[3];

    bool data_owned_ = false; //!< if yes, has to free the memory in data_map_
    datamap_t data_map_;

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
    /** @}*/

    const sid_t   level() const { return level_; }
    const real_t* hgrid() const { return hgrid_; }
    const real_t* xyz() const { return xyz_; }

    inline real_t xyz(const int id) const { return xyz_[id]; }
    inline real_t hgrid(const int id) const { return hgrid_[id]; }

    real_p       data(Field* fid);
    real_p       data(const Field* fid, const sid_t ida);
    const real_p data(const Field* fid) const;

    void AddField(Field* fid);
    void DeleteField(Field* fid);
    void AddFields(map<string, Field*>* fields) ;
};

/**
 * @brief pointer to an member function of the class @ref GridBlock
 */
using bop_t = void (GridBlock::*)(Field* fid);

#endif  // SRC_BLOCK_HPP_