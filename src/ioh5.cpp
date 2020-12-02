#include "ioh5.hpp"

#include "doop.hpp"
#include "mpi.h"

using std::numeric_limits;
using std::string;
using std::to_string;

/**
 * @brief add a string attribute to the loc hdf5 object
 * 
 * @note this function is inspired from AMReX
 * 
 * @param loc the object to add the attribute to
 * @param name the name of the attribute
 * @param data the string itself
 */
static void HDF5_AttributeString(hid_t loc, const char* name, const char* data) {
    hid_t  attr, atype, space;
    herr_t ret;

    space = H5Screate(H5S_SCALAR);
    atype = H5Tcopy(H5T_C_S1);
    H5Tset_size(atype, strlen(data) + 1);
    H5Tset_strpad(atype, H5T_STR_NULLTERM);
    attr = H5Acreate(loc, name, atype, space, H5P_DEFAULT, H5P_DEFAULT);
    m_assert(attr >= 0, "error while writting in hdf5");

    ret = H5Awrite(attr, atype, data);
    m_assert(ret >= 0, "error while writting in hdf5");

    H5Tclose(atype);
    H5Sclose(space);
    H5Aclose(attr);
}

/**
 * @brief Construct a new IOH5 given a folder. Any successive operator()() calls will dump into this folder
 * 
 * @param folder the folder where the I/O will be performed
 */
IOH5::IOH5(string folder) {
    m_begin;
    //-------------------------------------------------------------------------
    // store the folder
    folder_ = folder;
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief Destroy the IOH5 object
 * 
 */
IOH5::~IOH5() {
    m_begin;
    //-------------------------------------------------------------------------
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief dump the field with its name and an index = 0, see IOH5::operator()(ForestGrid* grid, const Field* field, const string name, const lid_t index)
 */
void IOH5::operator()(ForestGrid* grid, const Field* field) {
    m_begin;
    //-------------------------------------------------------------------------
    IOH5::operator()(grid, field, field->name(), 0);
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief dump the field with its name, see IOH5::operator()(ForestGrid* grid, const Field* field, const string name, const lid_t index)
 */
void IOH5::operator()(ForestGrid* grid, const Field* field, const lid_t index) {
    m_begin;
    //-------------------------------------------------------------------------
    IOH5::operator()(grid, field, field->name(), index);
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief dump the field with an index = 0, see IOH5::operator()(ForestGrid* grid, const Field* field, const string name, const lid_t index)
 */
void IOH5::operator()(ForestGrid* grid, const Field* field, const string name) {
    m_begin;
    //-------------------------------------------------------------------------
    IOH5::operator()(grid, field, name, 0);
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief dump the field (xdmf+h5) in the predefined folder, using the given name
 * 
 * A prefix `g_` is added ot the filename if the ghosts are dumped as well
 * 
 * @warning the ghosts of the field must be up to date as we IO min (M_N+1)^3 points, so that paraview handles it nicely
 * 
 * @param grid the grid supporting the field
 * @param field the field to dump
 * @param name the filename to use
 * @param index the index to append to the file's name
 */
void IOH5::operator()(ForestGrid* grid, const Field* field, const string name, const lid_t index) {
    m_begin;
    //-------------------------------------------------------------------------
    m_assert(field->ghost_status(), "the field <%s> has outdated ghosts, please update them, even if dump ghost is false",field->name().c_str());
    m_log("dumping field %s to disk (ghost = %d)", name.c_str(), dump_ghost_);
    // get the field name
    rank_t mpirank = grid->mpirank();
    rank_t rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    m_assert(mpirank == rank, "The rank from the forest does not correspond to the world one, this is not expected: %d vs %d", mpirank, rank);

    //................................................
    // get the filename
    char index_char[256];
    sprintf(index_char, "%6.6d", index);
    if (dump_ghost_) {
        filename_hdf5_ = "g_" + name + "_" + string(index_char) + ".h5";
        filename_xdmf_ = "g_" + name + "_" + string(index_char) + ".xmf";
    } else {
        filename_hdf5_ = name + "_" + string(index_char) + ".h5";
        filename_xdmf_ = name + "_" + string(index_char) + ".xmf";
    }

    // if the folder does not exist, create it
    struct stat st = {0};
    if (rank == 0 && stat(folder_.c_str(), &st) == -1) {
        mkdir(folder_.c_str(), 0770);  //create the folder if it does not exists
    }

    // get the bstride given the dumpghost status
    block_stride_ = (dump_ghost_) ? M_STRIDE : (M_N + 1);
    block_shift_  = (dump_ghost_) ? 0 : M_GS;

    // get the block offset and the number of global blocks
    size_t n_block_global = 0;
    size_t local_n_block  = (size_t)(grid->p4est_forest()->local_num_quadrants);
    MPI_Scan(&local_n_block, &block_offset_, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&local_n_block, &n_block_global, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    block_offset_ -= local_n_block;
    m_assert(sizeof(local_n_block) == sizeof(unsigned long), "the type used for the MPI Scan is not correct anymore");
    m_assert(sizeof(block_offset_) == sizeof(unsigned long), "the type used for the MPI Scan is not correct anymore");
    m_assert(sizeof(n_block_global) == sizeof(unsigned long), "the type used for the MPI Scan is not correct anymore");

    //................................................
    // print the header
    xmf_write_header_(grid, n_block_global, field->lda());
    hdf5_write_header_(grid, n_block_global, field->lda());

    //................................................
    // call first the xmf and then the h5 to reduce filesystem stress
    DoOpMesh(this, &IOH5::xmf_write_block_, grid, field);
    DoOpMesh(this, &IOH5::hdf5_write_block_, grid, field);

    //................................................
    // print the footer
    xmf_write_footer_(grid);
    hdf5_write_footer_(grid);
    // reset the dump_ghost to false for next dump
    dump_ghost_ = false;
    //-------------------------------------------------------------------------
    m_end;
}

// /**
//  * @brief Dumps the field contained in a given block
//  *
//  * @param qid the quadrant ID
//  * @param block the block to be dumped
//  * @param fid the field
//  */
// void IOH5::IOGridBlock(const qid_t* qid, GridBlock* block, const Field* fid) {
//     //-------------------------------------------------------------------------
//     xmf_write_block_(qid, block, fid);
//     hdf5_write_block_(qid, block, fid);
//     //-------------------------------------------------------------------------
// }

/**
 * @brief ask to dump the ghosts, this will be reset to false after the next dump
 * 
 * @param dump_ghost 
 */
void IOH5::dump_ghost(const bool dump_ghost) {
    dump_ghost_ = dump_ghost;
}

/**
 * @brief write the hdf5 header: create the file and initiate the datasets to contain the blocks (file and memory)
 * 
 * @param grid the grid
 * @param n_block_global the total number of block (global on the comm) to dump
 * @param lda the total dimension of the array to dump
 */
void IOH5::hdf5_write_header_(const ForestGrid* grid, const size_t n_block_global, const lda_t lda) {
    m_begin;
    //-------------------------------------------------------------------------
    //................................................
    // get the file opening properties
    hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(fapl, MPI_COMM_WORLD, MPI_INFO_NULL);
    // int alignment = block_stride_ * block_stride_ * block_stride_;
    // H5Pset_alignment(fapl, alignment, alignment);
    H5Pset_coll_metadata_write(fapl, true);
    H5Pset_all_coll_metadata_ops(fapl, true);
    hid_t dxpl = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_INDEPENDENT);

    // create the file ID
    const string filename = folder_ + string("/") + filename_hdf5_;
    hdf5_file_            = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
    m_assert(hdf5_file_ >= 0, "error while opening hdf5 file %s", filename.c_str());

    // close the properties
    H5Pclose(fapl);

    //................................................
    // write some usefull attribute data, just for fun
    const string filename_xmdf = folder_ + string("/") + filename_xdmf_;
    HDF5_AttributeString(hdf5_file_, "xmf file", filename_xmdf.c_str());
    HDF5_AttributeString(hdf5_file_, "written by ", "Murphy");

    //................................................
    // create the huuuge dataspace for all the blocks:
    hsize_t block_size[3] = {(hsize_t)(n_block_global * lda * block_stride_), (hsize_t)(block_stride_), (hsize_t)(block_stride_)};
    hdf5_dataspace_       = H5Screate_simple(3, block_size, NULL);
    m_assert(hdf5_dataspace_ >= 0, "error creating the dataspace");
    hdf5_dataset_ = H5Dcreate(hdf5_file_, "blocks", H5T_NATIVE_FLOAT, hdf5_dataspace_, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    m_assert(hdf5_dataset_ >= 0, "error creating the dataset");

    //................................................
    // create the memory layout: outer dim first!
    hsize_t memsize[3] = {(hsize_t)(M_STRIDE * lda), M_STRIDE, M_STRIDE};  // m_blockmemsize(lda);
    hdf5_memspace_     = H5Screate_simple(3, memsize, NULL);
    m_assert(hdf5_memspace_ >= 0, "error while creating the memory space");

    // get the needed hyperslab
    const hsize_t memstride[3] = {(hsize_t)block_stride_, (hsize_t)block_stride_, (hsize_t)block_stride_};  // distance between two blocks
    const hsize_t memblock[3]  = {(hsize_t)block_stride_, (hsize_t)block_stride_, (hsize_t)block_stride_};  // size of 1 block
    const hsize_t memcount[3]  = {1, 1, 1};                                                                 // number of blocks
    const hsize_t memoffset[3] = {(hsize_t)block_shift_, (hsize_t)block_shift_, (hsize_t)block_shift_};     // number of blocks

    // get the hyperslab
    herr_t status = H5Sselect_hyperslab(hdf5_memspace_, H5S_SELECT_SET, memoffset, memstride, memcount, memblock);
    m_assert(status >= 0, "failed to select the memory hyperslab");

    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief write down the block do the dataset
 * 
 * @warning the indexes, as used by hdf5, idx[3] are such that the first index (idx[0]) is the slowest rotating index
 * and the last one (idx[2]) is the fastest rotating index
 * 
 * @param qid 
 * @param block 
 * @param fid 
 */
void IOH5::hdf5_write_block_(const qid_t* qid, GridBlock* block, const Field* fid) {
    //-------------------------------------------------------------------------
    herr_t        status;  // error code
    const hsize_t field_lda = fid->lda();

    //................................................
    // set the constant file hyperslab params -> first index = slower rotating index
    const hsize_t filecount[3]  = {1, 1, 1};                                                                 // number of block
    const hsize_t filestride[3] = {(hsize_t)block_stride_, (hsize_t)block_stride_, (hsize_t)block_stride_};  // stride between two blocks
    const hsize_t fileblock[3]  = {(hsize_t)block_stride_, (hsize_t)block_stride_, (hsize_t)block_stride_};  // dimension of one block

    //looping over the vector components.
    for (lda_t ida = 0; ida < fid->lda(); ++ida) {
        // update the file offset
        const hsize_t fileoffset[3] = {((block_offset_ + qid->cid) * field_lda + ida) * block_stride_, 0, 0};

        // get the hyperslab
        status = H5Sselect_hyperslab(hdf5_dataspace_, H5S_SELECT_SET, fileoffset, filestride, filecount, fileblock);
        m_assert(status >= 0, "error while creating the file hyperslab");

        // do the writting for the data, we need to shift it by hand since memoffset cannot be < 0
        real_p data = block->pointer(fid, ida);
        status      = H5Dwrite(hdf5_dataset_, H5T_NATIVE_DOUBLE, hdf5_memspace_, hdf5_dataspace_, H5P_DEFAULT, data);
        m_assert(status >= 0, "error while creating the file hyperslab");
    }

    //-------------------------------------------------------------------------
}
/**
 * @brief write the footer of hdf5: close the spaces and the file
 * 
 * @param grid 
 */
void IOH5::hdf5_write_footer_(const ForestGrid* grid) {
    m_begin;
    //-------------------------------------------------------------------------
    H5Sclose(hdf5_memspace_);
    H5Sclose(hdf5_dataspace_);
    H5Dclose(hdf5_dataset_);
    H5Fclose(hdf5_file_);
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief write the header of the xmf + compute the offset of every rank in the file
 * 
 * @param grid the forest grid to be IOed
 * @param n_block_global the global number of block (total on the comm) involved in the IO
 * @param lda the number of dimensios to IO
 */
void IOH5::xmf_write_header_(const ForestGrid* grid, const size_t n_block_global, const lda_t lda) {
    m_begin;
    //-------------------------------------------------------------------------
    rank_t rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // check that the file name doesn't already exists
    struct stat st       = {0};
    string      filename = folder_ + string("/") + filename_xdmf_;
    m_assert(stat(filename.c_str(), &st) != 0, "ERROR while opening  <%s>, the file already exists", filename.c_str());
    // fopen the xmf, every proc
    int err = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_WRONLY | MPI_MODE_CREATE | MPI_MODE_EXCL, MPI_INFO_NULL, &xmf_file_);
    // if something went wrong, check if the file already exist of something else was baaad
    m_assert(err == MPI_SUCCESS, "ERROR while opening  <%s>, MPI_File_open failed (error = %d)", filename.c_str(), err);

    // the current position of current proc
    MPI_Offset offset = 0;
    MPI_File_seek(xmf_file_, offset, MPI_SEEK_SET);
    // write the header only in proc 0
    int header_count = 0;
    if (rank == 0) {
        // build the header
        const string header = string("<?xml version=\"1.0\" ?>\n") +
                              string("<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n") +
                              string("<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"2.2\">\n") +
                              string("<!-- Kindly generated by Murphy -->\n") +
                              string("  <Domain>\n") +
                              string("      <Grid GridType=\"Collection\" CollectionType=\"Spatial\">\n") +
                              string("<!-- put here your MR blocks -->");
        // get the len of the header (+1 for the \0 caracter)
        m_assert(header.size() < numeric_limits<int>::max(), "the header is too big for an int");
        header_count = header.size();

        // write the header
        MPI_Status status;
        MPI_File_write(xmf_file_, header.c_str(), header_count, MPI_CHAR, &status);
        // count the equivalent number of lines
        MPI_File_get_position(xmf_file_, &offset);
    }
    // do an empty quad string to check the length
    char msg[4096];
    memset(msg, 0, 4096);
    real_t zero[3] = {0.0, 0.0, 0.0};
    len_per_quad_  = xmf_core_(filename_hdf5_, zero, zero, 1, 1, 1, 1, 1, lda, 0, 0, 1, msg);

    // need to compute the shift of everybody
    size_t quad_len = grid->p4est_forest()->local_num_quadrants * len_per_quad_;
    size_t pos_end  = header_count + quad_len;

    // gt count, the position of proc i (in lines!!) as the sum of the position of procs 0 -> (i-1)
    size_t pos_curr = 0;
    MPI_Scan(&pos_end, &pos_curr, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    m_assert(sizeof(size_t) == sizeof(unsigned long), "the datatype is not correct!");
    // update the offset position
    offset = (pos_curr - quad_len);
    // shift the pointer of everybody to the write place
    MPI_File_seek(xmf_file_, offset, MPI_SEEK_SET);

    // compute the footer offset
    footer_offset_ = n_block_global * len_per_quad_ + header_count;

    // compute the length of the hdf5 file
    stride_global_ = n_block_global * lda * block_stride_;

#ifndef NDEBUG
    size_t len;
    MPI_Allreduce(&pos_end, &len, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    m_assert(sizeof(len) == sizeof(unsigned long), "this must be true");
    m_assert(sizeof(pos_end) == sizeof(unsigned long), "this must be true");
    if (rank == 0) {
        m_assert(len == (footer_offset_), "the length do not match: %ld vs %lld", len, footer_offset_);
    }
#endif

    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief write the xdmf footer using the precomputed @ref footer_offset_
 * 
 * @param grid 
 * @param n_block_global 
 */
void IOH5::xmf_write_footer_(const ForestGrid* grid) {
    m_begin;
    //-------------------------------------------------------------------------
    // we have to wait that everybody has finished before writting the footer
    rank_t rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
        // go to the end of the file
        MPI_File_seek(xmf_file_, footer_offset_, MPI_SEEK_SET);
        // write everything
        const string footer = string("\n<!-- end of your MR blocks -->\n") +
                              string("      </Grid>\n") +
                              string("  </Domain>\n") +
                              string("</Xdmf>\n");

        const size_t len = footer.size();

        // write the header
        MPI_Status status;
        MPI_File_write(xmf_file_, footer.c_str(), len, MPI_CHAR, &status);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    // finish!! close it
    MPI_File_close(&xmf_file_);
    //-------------------------------------------------------------------------
    m_end;
}
void IOH5::xmf_write_block_(const qid_t* qid, GridBlock* block, const Field* fid) {
    //-------------------------------------------------------------------------
    rank_t rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //#pragma omp critical
    char msg[4096];
    memset(msg, 0, 4096);
    size_t offset = (block_offset_ + qid->cid) * fid->lda() * block_stride_;
    size_t len    = xmf_core_(filename_hdf5_, block->hgrid(), block->xyz(), qid->tid, qid->qid, rank, block_stride_, M_GS - block_shift_, fid->lda(), offset, stride_global_, block->level(), msg);
    m_assert(len == len_per_quad_, "the len has changed, hence the file will be corrupted: now %ld vs stored %ld", len, len_per_quad_);
    // write the header
    MPI_Status status;
    MPI_File_write(xmf_file_, msg, len_per_quad_, MPI_CHAR, &status);
    //
}

size_t IOH5::xmf_core_(const string fname_h5, const real_t* hgrid, const real_t* xyz, const p4est_topidx_t tid, const p4est_locidx_t qid, const rank_t rank, const lid_t stride, const lid_t n_gs, const lda_t lda, const hsize_t offset, const hsize_t stride_global, const level_t level, char* msg) {
    //-------------------------------------------------------------------------
    // we need an extra space for the final character
    char line[256];
    memset(line, 0, 256);
    // - L1
    sprintf(line, "\n<!-- tree num %10.10d, quad num %10.10d for %10.10d -->\n", tid, qid, rank);
    strcat(msg, line);
    // - L2
    sprintf(line, "          <Grid Name=\"t%10.10d_r%10.10d_q%10.10d\" GridType=\"Uniform\">\n", tid, rank, qid);
    strcat(msg, line);
    // - L3
    sprintf(line, "              <Topology TopologyType=\"3DCoRectMesh\" Dimensions=\"%10.10d %10.10d %10.10d\"/>\n", stride, stride, stride);
    strcat(msg, line);
    // - L4
    sprintf(line, "              <Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n");
    strcat(msg, line);
    // - L5
    sprintf(line, "                  <DataItem Name=\"Origin\" Dimensions=\"3\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">\n");
    strcat(msg, line);
    // - L6
    sprintf(line, "                      %+10.8f %+10.8f %+10.8f\n", xyz[2] - n_gs * hgrid[2], xyz[1] - n_gs * hgrid[1], xyz[0] - n_gs * hgrid[0]);
    strcat(msg, line);
    // - L7
    sprintf(line, "                  </DataItem>\n");
    strcat(msg, line);
    // - L8
    sprintf(line, "                  <DataItem Name=\"Spacing\" Dimensions=\"3\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">\n");
    strcat(msg, line);
    // - L9
    sprintf(line, "                      %+10.8f %+10.8f %+10.8f\n", hgrid[2], hgrid[1], hgrid[0]);
    strcat(msg, line);
    // - L10
    sprintf(line, "                  </DataItem>\n");
    strcat(msg, line);
    // - L11
    sprintf(line, "              </Geometry>\n");
    strcat(msg, line);
    // - L12
    for (lda_t ida = 0; ida < lda; ++ida) {
        sprintf(line, "              <Attribute Name=\"data_%d\" AttributeType=\"Scalar\" Center=\"Node\">\n", ida);
        strcat(msg, line);
        // - L13
        sprintf(line, "                  <DataItem ItemType=\"HyperSlab\" Dimensions=\"%10.10d %10.10d %10.10d\" Type=\"HyperSlab\">\n", stride, stride, stride);
        strcat(msg, line);
        // - L13
        sprintf(line, "                     <DataItem Dimensions=\"3 3\" Format=\"XML\">\n");
        strcat(msg, line);
        // - L14
        sprintf(line, "                         %20.20lld %10.10d %10.10d\n", offset + ida * stride, 0, 0);
        strcat(msg, line);
        sprintf(line, "                         %10.10d %10.10d %10.10d\n", 1, 1, 1);
        strcat(msg, line);
        sprintf(line, "                         %10.10d %10.10d %10.10d\n", stride, stride, stride);
        strcat(msg, line);
        // - L15
        sprintf(line, "                     </DataItem>\n");
        strcat(msg, line);
        // - L13
        sprintf(line, "                     <DataItem Dimensions=\"%20.20lld %10.10d %10.10d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", stride_global, stride, stride);
        strcat(msg, line);
        // - L14
        sprintf(line, "                         %s:/blocks\n", fname_h5.c_str());
        strcat(msg, line);
        // - L15
        sprintf(line, "                     </DataItem>\n");
        strcat(msg, line);
        // - L15
        sprintf(line, "                  </DataItem>\n");
        strcat(msg, line);
        // - L16
        sprintf(line, "              </Attribute>\n");
        strcat(msg, line);
    }
#ifndef NDEBUG
    // - L17
    sprintf(line, "              <Attribute Name=\"rank\" AttributeType=\"Scalar\" Center=\"Grid\">\n");
    strcat(msg, line);
    // - L18
    sprintf(line, "                  <DataItem Dimensions=\"1\" NumberType=\"UInt\" Format=\"XML\">\n");
    strcat(msg, line);
    // - L19
    sprintf(line, "                      %10.10d\n", rank);
    strcat(msg, line);
    // - L20
    sprintf(line, "                  </DataItem>\n");
    strcat(msg, line);
    // - L21
    sprintf(line, "              </Attribute>\n");
    strcat(msg, line);
    // - L22
    sprintf(line, "              <Attribute Name=\"tree\" AttributeType=\"Scalar\" Center=\"Grid\">\n");
    strcat(msg, line);
    // - L23
    sprintf(line, "                  <DataItem Dimensions=\"1\" NumberType=\"UInt\" Format=\"XML\">\n");
    strcat(msg, line);
    // - L24
    sprintf(line, "                      %10.10d\n", tid);
    strcat(msg, line);
    // - L25
    sprintf(line, "                  </DataItem>\n");
    strcat(msg, line);
    // - L26
    sprintf(line, "              </Attribute>\n");
    strcat(msg, line);
#endif
    sprintf(line, "              <Attribute Name=\"level\" AttributeType=\"Scalar\" Center=\"Grid\">\n");
    strcat(msg, line);
    // - L23
    sprintf(line, "                  <DataItem Dimensions=\"1\" NumberType=\"UInt\" Format=\"XML\">\n");
    strcat(msg, line);
    // - L24
    sprintf(line, "                      %10.10d\n", level);
    strcat(msg, line);
    // - L25
    sprintf(line, "                  </DataItem>\n");
    strcat(msg, line);
    // - L26
    sprintf(line, "              </Attribute>\n");
    strcat(msg, line);
    // - L27
    sprintf(line, "         </Grid>");
    strcat(msg, line);
    //-------------------------------------------------------------------------
    return strlen(msg);
}