#include "ioh5.hpp"

#include "mpi.h"

IOH5::IOH5(string folder) {
    m_begin;
    //-------------------------------------------------------------------------
    // store the folder
    folder_ = folder;
    // create the line datatype required by MPI
    MPI_Type_contiguous(M_IOH5_LINE_LEN, MPI_CHAR, &line_type_);
    MPI_Type_commit(&line_type_);
    //-------------------------------------------------------------------------
    m_end;
}

IOH5::~IOH5() {
    m_begin;
    //-------------------------------------------------------------------------
    // free the created data_type
    MPI_Type_free(&line_type_);
    //-------------------------------------------------------------------------
    m_end;
}

void IOH5::DoOp(Grid* grid, Field* field) {
    m_begin;
    //-------------------------------------------------------------------------
    // get the field name
    filename_ = field->name();
    mpirank_  = grid->mpirank();
    // print the header
    xmf_write_header_(grid);
    hdf5_write_header_(grid);
    // // compute the ghost if needed:
    // if(dump_ghost_){
    //     exchange_t* ex= grid_ghost_start(grid,field);
    //     grid_ghost_end(grid,ex,field);
    // }
    // call the standard operator
    ConstOperatorF::DoOp(grid, field);
    // print the footer
    xmf_write_footer_(grid);
    hdf5_write_footer_(grid);
    //-------------------------------------------------------------------------
    m_end;
}

void IOH5::apply(const qid_t* qid, Block* block, const Field* fid) {
    m_begin;
    //-------------------------------------------------------------------------
    xmf_write_block_(qid, block, fid);
    hdf5_write_block_(qid, block, fid);
    //-------------------------------------------------------------------------
    m_end;
}

void IOH5::hdf5_write_header_(const Grid* grid) {
    m_begin;
    //-------------------------------------------------------------------------
    // get the full name
    string folder      = folder_ + "/" + filename_;
    string extFilename = folder_ + "/" + filename_ + "/" + filename_ + "_rank" + to_string(mpirank_) + ".h5";
    // create the folder if it does not exist
    struct stat st = {0};
    if (stat(folder.c_str(), &st) == -1) {
        mkdir(folder.c_str(), 0770);  //create the folder if it does not exists
    }
    // create the file ID
    hdf5_file_ = H5Fcreate(extFilename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    //-------------------------------------------------------------------------
    m_end;
}

void IOH5::hdf5_write_block_(const qid_t* qid, Block* block, const Field* fid) {
    m_begin;
    //-------------------------------------------------------------------------
    // get the attribute
    string attribute = "tree" + to_string(qid->tid) + "_quad" + to_string(qid->qid);
    hid_t  filespace;  // dataspaces
    hid_t  fileset;    // datasets
    hid_t  memspace;
    herr_t status;  // error code

    //-------------------------------------------------------------------------
    /** - Create the file dataspace and dataset  */
    /** \warning In the dataspace, the last index must be the index of the vector 
     * component (requirement from xdmf). However, in memory, the index of the
     * vector component is the first one. We will thus need to fill the file 
     * with the data, component by component, and using a stride of lda.
    //-----------------------------------------------------------------------*/
    // the file information is given by the global size = total size reserved for the file
    hsize_t field_dims[4] = {M_N, M_N, M_N, (hsize_t)fid->lda()};
    if (dump_ghost_) {
        field_dims[0] = M_STRIDE;
        field_dims[1] = M_STRIDE;
        field_dims[2] = M_STRIDE;
    }

    // setup the options to have a chuncked dataset
    hid_t   plist_id   = H5Pcreate(H5P_DATASET_CREATE);
    hsize_t chk_dim[4] = {4, 4, 4, 1};
    H5Pset_chunk(plist_id, 4, chk_dim);

    // create dataset and dataspace = the whole hard memory reserved for the file
    filespace = H5Screate_simple(4, field_dims, NULL);
    fileset   = H5Dcreate(hdf5_file_, attribute.c_str(), H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(filespace);
    H5Pclose(plist_id);

    //-------------------------------------------------------------------------
    /** - select the hyperslab inside the file dataset (=writting location)  */
    //-------------------------------------------------------------------------
    // set the full memory space
    hsize_t memsize[4] = {(hsize_t)fid->lda(), M_STRIDE, M_STRIDE, M_STRIDE};  // the full memory
    memspace           = H5Screate_simple(4, memsize, NULL);

    // caracteristics of the hyperslab in the file
    hsize_t count[4]  = {1, 1, 1, 1};                    // how many blocks to write
    hsize_t stride[4] = {1, 1, 1, (hsize_t)fid->lda()};  // distance between 2 blocks
    hsize_t bsize[4]  = {M_N, M_N, M_N, 1};              // the block size = the local size
    hsize_t offset[4] = {0, 0, 0, 0};
    // caracteristics fo the hyperslabs in the
    hsize_t memcount[4]  = {1, M_N, M_N, M_N};
    hsize_t memblock[4]  = {1, 1, 1, 1};
    hsize_t memstride[4] = {(hsize_t)fid->lda(), 1, 1, 1};
    hsize_t memoffset[4] = {0, M_GS, M_GS, M_GS};  // offset in memory

    if (dump_ghost_) {
        // size fo the file block
        bsize[0] = M_STRIDE;
        bsize[1] = M_STRIDE;
        bsize[2] = M_STRIDE;
        // number of memory blocks
        memcount[1] = M_STRIDE;
        memcount[2] = M_STRIDE;
        memcount[3] = M_STRIDE;
        // offset in memory
        memoffset[1] = 0;
        memoffset[2] = 0;
        memoffset[3] = 0;
    }

    //looping over the vector components.
    for (int lia = 0; lia < fid->lda(); lia++) {
        // set the offset to match the current dir
        offset[3]    = (hsize_t)lia;
        memoffset[0] = (hsize_t)lia;
        // select the hyperslab of the file
        filespace = H5Dget_space(fileset);
        status    = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, stride, count, bsize);
        // select the hyperslab of the memory
        status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, memoffset, memstride, memcount, memblock);
        // do the writting
        status = H5Dwrite(fileset, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, block->data(fid));
    }

    //-------------------------------------------------------------------------
    /** - close everything */
    //-------------------------------------------------------------------------
    H5Dclose(fileset);
    H5Sclose(filespace);
    H5Sclose(memspace);
    //-------------------------------------------------------------------------
    m_end;
}

void IOH5::hdf5_write_footer_(const Grid* grid) {
    H5Fclose(hdf5_file_);
}

void IOH5::xmf_write_header_(const Grid* grid) {
    m_begin;
    //-------------------------------------------------------------------------
    sc_MPI_Comm comm  = grid->mpicomm();
    string      fname = folder_ + "/" + filename_ + ".xmf";
    // rank 0 creates the folder if it doesn't exist
    if (mpirank_ == 0) {
        struct stat st = {0};
        if (stat(folder_.c_str(), &st) == -1) {
            mkdir(folder_.c_str(), 0770);
        }
    }
    // fopen the xmf, every proc
    int err = MPI_File_open(comm, fname.c_str(), MPI_MODE_WRONLY | MPI_MODE_CREATE | MPI_MODE_EXCL, MPI_INFO_NULL, &xmf_file_);
    if (err != MPI_SUCCESS) {
        m_log("ERROR while opening  <%s>, the file may be corrupted", fname.c_str());
    }
    // the current position of current proc
    MPI_Offset offset = 0;
    // go to the begining of the file
    MPI_File_seek(xmf_file_, offset, MPI_SEEK_SET);
    // write the header only in proc 0
    if (mpirank_ == 0) {
        // write the stuff and count
        char msg[7 * M_IOH5_LINE_LEN + 1];
        sprintf(msg, "%-128s", "<?xml version=\"1.0\" ?>");
        sprintf(msg, "%s%-128s", msg, "\n<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>");
        sprintf(msg, "%s%-128s", msg, "\n<Xdmf xmlns:xi=\"http://www.w3.org/2003/XInclude\" Version=\"3.0\">");
        sprintf(msg, "%s%-128s", msg, "\n<!-- Kindly generated by Murphy -->");
        sprintf(msg, "%s%-128s", msg, "\n  <Domain>");
        sprintf(msg, "%s%-128s", msg, "\n      <Grid GridType=\"Collection\" CollectionType=\"Spatial\">");
        sprintf(msg, "%s%-128s", msg, "\n<!-- put here your MR blocks -->");

        // write the header
        MPI_Status status;
        MPI_File_write(xmf_file_, msg, 7, line_type_, &status);
        // count the equivalent number of lines
        MPI_File_get_position(xmf_file_, &offset);
    }
    // the number of quadrant written by each proc (in bytes)
    int count  = offset / (M_IOH5_LINE_LEN * sizeof(char));  // number of lines in the header
    int qnum   = grid->local_num_quadrants() * 17;           // number of lines per quadrant
    int posend = count + qnum;                               // the final position
    // get count, the position of proc i (in lines!!) as the sum of the position of procs 0 -> (i-1)
    MPI_Scan(&posend, &count, 1, MPI_INT, MPI_SUM, comm);
    // update the offset position
    offset = (count - qnum) * (M_IOH5_LINE_LEN * sizeof(char));
    // shift the pointer of everybody to the write place
    MPI_File_seek(xmf_file_, offset, MPI_SEEK_SET);
    //-------------------------------------------------------------------------
    m_end;
}

void IOH5::xmf_write_footer_(const Grid* grid) {
    m_begin;
    //-------------------------------------------------------------------------
    // we have to wait that everybody has finished before writting the footer
    MPI_Barrier(grid->mpicomm());
    if (mpirank_ == 0) {
        // go to the end of the file
        MPI_File_seek(xmf_file_, 0, MPI_SEEK_END);
        // write everything
        char msg[512];
        sprintf(msg, "%s", "\n<!-- end of your MR blocks -->");
        sprintf(msg, "%s%s", msg, "\n      </Grid>");
        sprintf(msg, "%s%s", msg, "\n  </Domain>");
        sprintf(msg, "%s%s", msg, "\n</Xdmf>");

        // write the header
        MPI_Status status;
        size_t     msglen = strlen(msg);
        m_assert(msglen <= 512, "the provided length is not large enough");
        MPI_File_write(xmf_file_, msg, msglen, MPI_CHAR, &status);
    }
    // finish!! close it
    MPI_File_close(&xmf_file_);
    //-------------------------------------------------------------------------
    m_end;
}
void IOH5::xmf_write_block_(const qid_t* qid, Block* block, const Field* fid) {
    m_begin;
    //-------------------------------------------------------------------------
    string fname     = filename_ + "_";
    string attribute = "tree" + to_string(qid->tid) + "_quad" + to_string(qid->qid);
    string folder    = folder_;
    // get grid specs
    const real_t* hgrid = block->hgrid();
    const real_t* xyz   = block->xyz();

    // everyblock is a uniform grid of fixed size
    // the size of the grid has to be +1 since the 3DCoRectMesh is defined as vertex-centered
    char line[M_IOH5_LINE_LEN];
    // we need an extra space for the final character
    char msg[17 * M_IOH5_LINE_LEN + 1];
    // - L1
    sprintf(line, "\n<!-- tree num %d, quad num %d for %d -->", qid->tid, qid->qid, mpirank_);
    sprintf(msg, "%-128s", line);
    // - L2
    sprintf(line, "\n          <Grid Name=\"%s\" GridType=\"Uniform\">", attribute.c_str());
    sprintf(msg, "%s%-128s", msg, line);
    // - L3
    if (dump_ghost_) {
        sprintf(line, "\n              <Topology TopologyType=\"3DCoRectMesh\" Dimensions=\"%d %d %d\"/>", M_STRIDE + 1, M_STRIDE + 1, M_STRIDE + 1);
    } else {
        sprintf(line, "\n              <Topology TopologyType=\"3DCoRectMesh\" Dimensions=\"%d %d %d\"/>", M_N + 1, M_N + 1, M_N + 1);
    }
    sprintf(msg, "%s%-128s", msg, line);
    // - L4
    sprintf(msg, "%s%-128s", msg, "\n              <Geometry GeometryType=\"ORIGIN_DXDYDZ\">");
    // - L5
    sprintf(msg, "%s%-128s", msg, "\n                  <DataItem Name=\"Origin\" Dimensions=\"3\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">");
    // - L6
    if (dump_ghost_) {
        sprintf(line, "\n                      %10.8f %10.8f %10.8f", xyz[2] - M_GS * hgrid[2], xyz[1] - M_GS * hgrid[1], xyz[0] - M_GS * hgrid[0]);
    } else {
        sprintf(line, "\n                      %10.8f %10.8f %10.8f", xyz[2], xyz[1], xyz[0]);
    }
    sprintf(msg, "%s%-128s", msg, line);
    // - L7
    sprintf(msg, "%s%-128s", msg, "\n                  </DataItem>");
    // - L8
    sprintf(msg, "%s%-128s", msg, "\n                  <DataItem Name=\"Spacing\" Dimensions=\"3\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">");
    // - L9
    sprintf(line, "\n                      %10.8f %10.8f %10.8f", hgrid[2], hgrid[1], hgrid[0]);
    sprintf(msg, "%s%-128s", msg, line);
    // - L10
    sprintf(msg, "%s%-128s", msg, "\n                  </DataItem>");
    // - L11
    sprintf(msg, "%s%-128s", msg, "\n              </Geometry>");
    // - L12
    sprintf(line, "\n              <Attribute Name=\"data\" AttributeType=\"Vector\" Center=\"Cell\">");
    sprintf(msg, "%s%-128s", msg, line);
    // - L13
    if (dump_ghost_) {
        sprintf(line, "\n                  <DataItem Dimensions=\"%d %d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">", M_STRIDE, M_STRIDE, M_STRIDE, fid->lda());
    } else {
        sprintf(line, "\n                  <DataItem Dimensions=\"%d %d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">", M_N, M_N, M_N, fid->lda());
    }

    sprintf(msg, "%s%-128s", msg, line);
    // - L14
    sprintf(line, "\n                      %s/%s_rank%d.h5:/%s", filename_.c_str(), filename_.c_str(), mpirank_, attribute.c_str());
    sprintf(msg, "%s%-128s", msg, line);
    // - L15
    sprintf(msg, "%s%-128s", msg, "\n                  </DataItem>");
    // - L16
    sprintf(msg, "%s%-128s", msg, "\n              </Attribute>");
    // - L17
    sprintf(msg, "%s%-128s", msg, "\n         </Grid>");
    // write the header
    MPI_Status status;
    MPI_File_write(xmf_file_, msg, 17, line_type_, &status);
    //-------------------------------------------------------------------------
    m_end;
}