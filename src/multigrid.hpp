#ifndef SRC_MULTIGRID_HPP
#define SRC_MULTIGRID_HPP

#include <map>
#include <string>

#include "field.hpp"
#include "grid.hpp"
#include "mgfamily.hpp"
#include "murphy.hpp"
#include "p8est.h"
#include "partitioner.hpp"

class Multigrid {
   protected:
    sid_t eta_1_ = 3;
    sid_t eta_2_ = 3;
    real_t alpha_ = 1.95;

    sid_t fft_level_ = 0;
    sid_t max_level_ = 0;
    sid_t n_level_   = 0;

    Grid**        grids_;
    Partitioner** parts_;

    map<string, Field*> map_fields_;

    lid_t*     n_family_member_;
    MGFamily** families_;

    sid_t tmp_ilevel_ = 0;

   public:
    Multigrid(Grid* grid, sid_t fft_level, Field* src, Field* trg, Field* res);
    ~Multigrid();

    void Solve();

    inline sid_t                curr_level() const { return fft_level_ + tmp_ilevel_; }
    inline MGFamily*            curr_family() { return families_[tmp_ilevel_]; }
    inline map<string, Field*>* map_fields() { return &(map_fields_); }
};

#endif  // SRC_MULTIGRID_HPP
