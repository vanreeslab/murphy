#ifndef SRC_IOH5_HPP
#define SRC_IOH5_HPP

// c headers
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include <iostream>

#include "forestgrid.hpp"
#include "hdf5.h"
#include "murphy.hpp"
#include "operator.hpp"

#define M_IOH5_LINE_LEN 128

using std::string;

/**
 * @brief define the I/O opertion
 *
 * The current I/O relies on a common XMF file, written by everybody and a sequential hdf5 I/O, i.e. one file per proc.
 * This is NOT a long term solution and it has to change, hence a restricted documentation is provided
 * 
 */
class IOH5 : public ConstOperatorF {
   protected:
    int  mpirank_    = 0;
    bool dump_ghost_ = false;

    string folder_;
    string filename_;

    hid_t        hdf5_file_;
    MPI_File     xmf_file_;
    MPI_Datatype line_type_;

    // headers
    void xmf_write_header_(const ForestGrid* grid);
    void hdf5_write_header_(const ForestGrid* grid);
    // footers
    void xmf_write_footer_(const ForestGrid* grid);
    void hdf5_write_footer_(const ForestGrid* grid);
    // blocks
    void xmf_write_block_(const qid_t* qid, GridBlock* block, const Field* fid);
    void hdf5_write_block_(const qid_t* qid, GridBlock* block, const Field* fid);
    // write a block

   public:
    IOH5(string folder);
    ~IOH5();

    void dump_ghost(const bool dump_ghost);

    void ApplyConstOpF(const qid_t* qid, GridBlock* block, const Field* fid) override;

    void operator()(ForestGrid* grid, Field* field, string name);
    void operator()(ForestGrid* grid, Field* field) override;
};

#endif