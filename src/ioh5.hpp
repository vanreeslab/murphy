#ifndef SRC_IOH5_HPP
#define SRC_IOH5_HPP

// c headers
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include <iostream>

#include "field.hpp"
#include "forestgrid.hpp"
#include "gridblock.hpp"
#include "hdf5.h"
#include "murphy.hpp"

#define M_IOH5_LINE_LEN 128

/**
 * @brief perform the I/O using HDF5 and XDMF files
 * 
 * @note The vector fields will be represented as scalars.
 * 
 * The IO will appear in Paraview as a MultiBlock dataset. In order to enable most of the filters, you have to options:
 * 
 * - the `Merge Blocks` filter will create an `Unstructured` dataset. That enables almost all the filters (except the volume rending!!) while preserving the info quality.
 * However, it might require a large amount of memory to perfom!
 * 
 * - the `Resample to Image` filter will downsample the dataset to a uniform representation. This is currently the only way to use the volume rendering filter.
 * The amount of needed memory is smaller but so is the image quality.
 * 
 */
class IOH5 {
   protected:
    bool   dump_ghost_    = false;  //!< true if we io the ghost points with us
    lid_t  block_stride_  = 0;      //!< the stride of one block, depends on dump_ghost_
    lid_t  block_shift_   = 0;      //!< the shift in memory to apply on the block ptr, depends on dump_ghost_
    size_t block_offset_  = 0;      //!< the offset in the xdmf file for the current rank
    size_t stride_global_ = 0;      //!< total length of the hdf5 dataset, must be indicated in the xdmf of every block.

    std::string folder_;         //!< the folder where do dump the files
    std::string filename_hdf5_;  //!< the filename of the hdf5 file
    std::string filename_xdmf_;  //!< the filenanem of the xdmf file

    hid_t hdf5_file_;       //!< hdf5 object refering to the filename
    hid_t hdf5_dataset_;    //!< hdf5 dataset to store the blocks
    hid_t hdf5_dataspace_;  //!< hdf5 dataspace to store the blocks, is changed for every dimension of every block
    hid_t hdf5_memspace_;   //!< hdf5 dataspace for the memory, if fixed as one dimension in a bloc

    size_t     len_per_quad_;   //!< the xdmf string length per bloc, used to compute the offest and re
    MPI_File   xmf_file_;       //!< xdmf distributed file
    MPI_Offset footer_offset_;  //<! offset of the footer, known by everybody but used only by rank 0.

    // headers
    void hdf5_write_header_(const ForestGrid* grid, const size_t n_block_global, const lda_t lda);
    void xmf_write_header_(const ForestGrid* grid, const size_t n_block_global, const lda_t lda);

    // footers
    void xmf_write_footer_(const ForestGrid* grid);
    void hdf5_write_footer_(const ForestGrid* grid);
    // blocks
    void xmf_write_block_(const qid_t* qid, GridBlock* block, const Field* fid);
    void hdf5_write_block_(const qid_t* qid, GridBlock* block, const Field* fid);
    // write a block

   public:
    explicit IOH5(std::string folder);
    ~IOH5();

    void dump_ghost(const bool dump_ghost);

    void operator()(ForestGrid* grid, const Field* field, const std::string name);
    void operator()(ForestGrid* grid, const Field* field);

   protected:
    size_t xmf_core_(const std::string fname_h5, const real_t* hgrid, const real_t* xyz, const p4est_topidx_t tid, const p4est_locidx_t qid, const rank_t rank, const lid_t stride, const lid_t n_gs, const lda_t lda, const hsize_t offset, const hsize_t stride_global, const level_t level, char* msg);
};

#endif