#ifndef SRC_IOHDF5_HPP
#define SRC_IOHDF5_HPP

// c headers
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include <iostream>

#include "forestgrid.hpp"
#include "hdf5.h"
#include "murphy.hpp"
#include "operator.hpp"

class IOHDF5 : public ConstOperatorF {
   protected:
    bool   dump_ghost_ = false;
    string folder_;
    string filename_;

    void hdf5_write_header_(const ForestGrid* grid, const Field* fid, const level_t min_lvl, const level_t max_lvl);
    void hdf5_write_block_(const qid_t* qid, GridBlock* block, const Field* fid);

   public:
    IOHDF5(string folder);
    ~IOHDF5();

    void dump_ghost(const bool dump_ghost);

    void ApplyConstOpF(const qid_t* qid, GridBlock* block, const Field* fid) override;
    void operator()(ForestGrid* grid, Field* field, std::string name);
    void operator()(ForestGrid* grid, Field* field) override;

   protected:
    void HDF5_WriteMetaData_(hid_t fid, const Field* field, const level_t min_lvl, const level_t max_lvl, const real_t L[3]);
};

#endif  // SRC_IOHDF5_HPP