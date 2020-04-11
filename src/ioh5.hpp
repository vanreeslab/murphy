#ifndef SRC_IOH5_HPP
#define SRC_IOH5_HPP

// c headers
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include <iostream>

#include "hdf5.h"
#include "murphy.hpp"
#include "operator.hpp"

// define macros to strigyfy, both are required!
#define STR(a) ZSTR(a)
#define ZSTR(a) #a

#define M_IOH5_LINE_LEN 128

using std::string;
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
    void xmf_write_header_(const Grid* grid);
    void hdf5_write_header_(const Grid* grid);
    // footers
    void xmf_write_footer_(const Grid* grid);
    void hdf5_write_footer_(const Grid* grid);
    // blocks
    void xmf_write_block_(const qid_t* qid, GridBlock* block, const Field* fid);
    void hdf5_write_block_(const qid_t* qid, GridBlock* block, const Field* fid);
    // write a block

   public:
    IOH5(string folder);
    ~IOH5();

    void dump_ghost(const bool dump_ghost) { dump_ghost_ = dump_ghost; };

    void apply(const qid_t* qid, GridBlock* block, const Field* fid) override;
    void operator()(Grid* grid, Field* field) override;
};

#endif