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

/**
 * @brief define the I/O using HDF5 and XDMF
 * 
 */
class IOH5 : public ConstOperatorF {
   protected:
    bool    dump_ghost_     = false;  //!< true if we io the ghost points with us
    hsize_t block_stride_   = 0;
    hsize_t block_shift_    = 0;
    hsize_t block_offset_   = 0;
    hsize_t n_block_global_ = 0;

    string folder_;
    string filename_hdf5_;
    string filename_xdmf_;

    hid_t hdf5_file_;
    hid_t hdf5_dataset_;
    hid_t hdf5_dataspace_;
    hid_t hdf5_memspace_;

    MPI_File xmf_file_;
    size_t   len_per_quad_;
    size_t   footer_offset_;

    // headers
    void hdf5_write_header_(const ForestGrid* grid, const hsize_t n_block_global, const lda_t lda);
    void xmf_write_header_(const ForestGrid* grid, const hsize_t n_block_global, const lda_t lda);

    // footers
    void xmf_write_footer_(const ForestGrid* grid, const size_t n_block_global);
    void hdf5_write_footer_(const ForestGrid* grid);
    // blocks
    void xmf_write_block_(const qid_t* qid, GridBlock* block, const Field* fid);
    void hdf5_write_block_(const qid_t* qid, GridBlock* block, const Field* fid);
    // write a block

   public:
    explicit IOH5(string folder);
    ~IOH5();

    void dump_ghost(const bool dump_ghost);

    void ApplyConstOpF(const qid_t* qid, GridBlock* block, const Field* fid) override;

    void operator()(ForestGrid* grid, Field* field, string name);
    void operator()(ForestGrid* grid, Field* field) override;

    protected:
    size_t xmf_core_(const string fname_h5, const real_t* hgrid, const real_t* xyz, const p4est_topidx_t tid, const p4est_locidx_t qid, const rank_t rank, const lid_t stride, const lid_t n_gs, const lda_t lda, const hsize_t offset, const hsize_t stride_global, char* msg);
};

#endif