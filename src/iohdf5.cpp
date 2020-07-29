#include "iohdf5.hpp"
#include "mpi.h"


#include <cstdarg>

using std::to_string;
using std::string;
using ullong = unsigned long long;

// inspired from AMRex - file AMReX_PlotFileUtil.cpp
static void HDF5_AttributeReal(hid_t loc, const char *name, hsize_t n, const real_t *data)
{
    herr_t ret;
    hid_t   attr, attr_space;
    hsize_t dims = n;
    // create the attribute space
    attr_space = H5Screate_simple(1, &dims, NULL);
    attr       = H5Acreate(loc, name, H5T_NATIVE_DOUBLE, attr_space, H5P_DEFAULT, H5P_DEFAULT);
    m_assert(attr >= 0, "error while writting in hdf5");
    // write it
    ret = H5Awrite(attr, H5T_NATIVE_DOUBLE, (void*)data);
    m_assert(ret >= 0, "error while writting in hdf5");
    H5Sclose(attr_space);
    H5Aclose(attr);
}
// inspired from AMRex
static void HDF5_AttributeInt(hid_t loc, const char* name, hsize_t n, const int* data) {
    herr_t  ret;
    hid_t   attr, attr_space;
    hsize_t dims = n;
    // create the attribute space
    attr_space = H5Screate_simple(1, &dims, NULL);
    attr       = H5Acreate(loc, name, H5T_NATIVE_INT, attr_space, H5P_DEFAULT, H5P_DEFAULT);
    m_assert(attr >= 0, "error while writting in hdf5");
    // write it
    ret = H5Awrite(attr, H5T_NATIVE_INT, (void*)data);
    m_assert(ret >= 0, "error while writting in hdf5");
    H5Sclose(attr_space);
    H5Aclose(attr);
}
// inspired from AMRex
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

static int HDF5_DatasetReal(hid_t loc, const char* name, hsize_t n, const real_t* data) {
    herr_t  ret;
    hid_t   dset, dset_space;
    hsize_t dims = n;

    dset_space = H5Screate_simple(1, &dims, NULL);

    dset = H5Dcreate(loc, name, H5T_NATIVE_FLOAT, dset_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    m_assert(dset >= 0, "error while writting in HDF5");

    ret = H5Dwrite(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, (void*)data);
    m_assert(ret >= 0, "error while writting in HDF5");
    H5Sclose(dset_space);
    H5Aclose(dset);
    return 1;
}

/**
 * @brief Construct a new IOH5 given a folder. Any successive operator()() calls will dump into this folder
 * 
 * @param folder 
 */
IOHDF5::IOHDF5(string folder) {
    m_begin;
    //-------------------------------------------------------------------------
    // store the folder
    folder_ = folder;
    //-------------------------------------------------------------------------
    m_end;
}

IOHDF5::~IOHDF5() {
    m_begin;
    //-------------------------------------------------------------------------
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief dump the field in the predefined folder, using the name fo the field for the filename
 * 
 * @param grid 
 * @param field 
 */
void IOHDF5::operator()(ForestGrid* grid, Field* field) {
    m_begin;
    //-------------------------------------------------------------------------
    IOHDF5::operator()(grid, field, field->name());
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief dump the field in the predefined folder, using the given name for the filename and reset the dump_ghost_ variable to false
 * 
 * @param grid the grid on which the field lives
 * @param field the field to dump
 * @param name the filename to use (`name.xmf`, `name/name_ranki.hdf5`)
 */
void IOHDF5::operator()(ForestGrid* grid, Field* field, string name) {
    m_begin;
    //-------------------------------------------------------------------------
    // get the file name
    if (dump_ghost_) {
        filename_ = folder_ + "/g_" + name + ".h5";
    } else {
        filename_ = folder_ + "/" + name + ".h5";
    }
    // create the folder if it does not exist
    rank_t rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        struct stat st = {0};
        if (stat(folder_.c_str(), &st) == -1) {
            mkdir(folder_.c_str(), 0770);  //create the folder if it does not exists
        }
    }
    m_log("dumping field %s to disk (ghost = %d)", name.c_str(), dump_ghost_);

    // get the forest and mesh info
    p8est_mesh_t* mesh   = grid->mesh();
    p8est_t*      forest = grid->forest();

    // given the dump_ghost status, get how much is 1 block:
    const hsize_t block_stride = ((dump_ghost_) ? (M_STRIDE) : (M_N));
    const hsize_t block_mem    = block_stride * block_stride * block_stride;
    const hsize_t field_lda    = (hsize_t)field->lda();

    //................................................
    // create the file and write the metadata
    level_t  max_lvl, min_lvl;
    level_t  local_max_lvl = -1, local_min_lvl = P8EST_MAXLEVEL;
    iblock_t n_block_old_level = 0;
    for (level_t il = 0; il < P8EST_MAXLEVEL; ++il) {
        iblock_t n_block_curr_level = p4est_NumQuadOnLevel(mesh, il);
        // increment the max if we have ghost and the min if we had nothing and we now have blocks
        local_max_lvl     = (n_block_curr_level > 0) ? il : local_max_lvl;
        local_min_lvl     = (n_block_curr_level > 0 && n_block_old_level == 0) ? il : local_min_lvl;
        n_block_old_level = n_block_curr_level;
        m_log("I have %d block on level %d -> %d to %d", n_block_curr_level, il, local_min_lvl, local_max_lvl);
    }
    // local_min_lvl++;  // the min_lvl = the first non-zero
    m_assert(sizeof(level_t) == 1, "please change the reduce datatype!");
    MPI_Allreduce(&local_min_lvl, &min_lvl, 1, MPI_INT8_T, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&local_max_lvl, &max_lvl, 1, MPI_INT8_T, MPI_MAX, MPI_COMM_WORLD);
    hdf5_write_header_(grid, field, min_lvl, max_lvl);
    m_log("levels going from %d to %d", min_lvl, max_lvl);

    //................................................
    // create the needed datatype for each level
    hid_t center_id = H5Tcreate(H5T_COMPOUND, 3 * sizeof(int));
    hid_t babox_id  = H5Tcreate(H5T_COMPOUND, 6 * sizeof(int));
    H5Tinsert(babox_id, "lo_i", 0 * sizeof(int), H5T_NATIVE_INT);
    H5Tinsert(babox_id, "lo_j", 1 * sizeof(int), H5T_NATIVE_INT);
    H5Tinsert(babox_id, "lo_k", 2 * sizeof(int), H5T_NATIVE_INT);
    H5Tinsert(babox_id, "hi_i", 3 * sizeof(int), H5T_NATIVE_INT);
    H5Tinsert(babox_id, "hi_j", 4 * sizeof(int), H5T_NATIVE_INT);
    H5Tinsert(babox_id, "hi_k", 5 * sizeof(int), H5T_NATIVE_INT);
    H5Tinsert(center_id, "i", 0 * sizeof(int), H5T_NATIVE_INT);
    H5Tinsert(center_id, "j", 1 * sizeof(int), H5T_NATIVE_INT);
    H5Tinsert(center_id, "k", 2 * sizeof(int), H5T_NATIVE_INT);

    // opening the file
    hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(fapl, MPI_COMM_WORLD, MPI_INFO_NULL);
    int alignment = 16 * 1024 * 1024;
    H5Pset_alignment(fapl, alignment, alignment);
    H5Pset_coll_metadata_write(fapl, true);
    H5Pset_all_coll_metadata_ops(fapl, true);

    hid_t dxpl = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_INDEPENDENT);

    // All process open the file
    hid_t hdf5_file = H5Fopen(filename_.c_str(), H5F_ACC_RDWR, fapl);
    m_assert(hdf5_file >= 0, "error while opening hdf5 file %s", filename_.c_str());

    //................................................
    // create the needed block dataspace as a single 1D vector, which is always the same
    // memory layout, no space needed, it's the memory layout, define one line as the mem is not always continuous
    hsize_t memsize      = m_blockmemsize(field_lda);
    hid_t   memdataspace = H5Screate_simple(1, &memsize, NULL);
    // set the selected space to empty
    herr_t status;
    status = H5Sselect_none(memdataspace);
    m_assert(status >= 0, "failed to empty the dataspace");
    // get the needed memory layout
    const hsize_t memstride = M_STRIDE;      // distance between two blocks
    const hsize_t memblock  = block_stride;  // size of 1 block
    const hsize_t memcount  = block_stride;  // number of blocks
    for (lda_t ida = 0; ida < field_lda; ++ida) {
        for (lid_t id2 = 0; id2 < block_stride; id2++) {
            // this is the offset in Z, computed based on the full stride, not the block one
            const hsize_t memoffset = memstride * memstride * (id2 + memstride * ida);
            // select the correct hyperslab
            status = H5Sselect_hyperslab(memdataspace, H5S_SELECT_OR, &memoffset, &memstride, &memcount, &memblock);
            m_assert(status >= 0, "failed to select the memory hyperslab");
        }
    }
    // store the shift in memory
    const size_t local_mem_offset = (dump_ghost_) ? 0 : m_szeroidx(0, M_GS, M_STRIDE);

    //................................................
    // every cpu has to go through every level, not only the local ones because all the metadata MUST be globally created
    for (level_t il = min_lvl; il <= max_lvl; ++il) {
        // everybody opens the level group
        char level_name[32];
        sprintf(level_name, "level_%d", (il - min_lvl));
        hid_t grp = H5Gopen(hdf5_file, level_name, H5P_DEFAULT);
        m_assert(grp >= 0, "fail to open level group %d", il);

        // get how many local blocks are involved
        iblock_t n_block_local = p4est_NumQuadOnLevel(mesh, il);
        // compute the block offset on each rank and how many blocks on the current level in total
        iblock_t n_block_offset, n_block_global;
        MPI_Allreduce(&n_block_local, &n_block_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Scan(&n_block_local, &n_block_offset, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        // remember if I am the last cpu + small correction on the block offset
        const bool is_last = (n_block_offset == n_block_global) && (n_block_global > 0);
        n_block_offset -= n_block_local;

        // complete the metadatas
        int ngrid = 1;
        HDF5_AttributeInt(grp, "ngrid", 1, &ngrid);

        //................................................
        // the datas, setup to the global size of the level
        string        dataname("data:datatype=0");
        const hsize_t level_size = n_block_global * block_mem * field_lda;
        hid_t         dataspace  = H5Screate_simple(1, &level_size, NULL);
        hid_t         dataset    = H5Dcreate(grp, dataname.c_str(), H5T_NATIVE_FLOAT, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        const hsize_t filestride = block_mem * field_lda;  // stride between two blocks
        const hsize_t fileblock  = block_mem * field_lda;  // block size
        const hsize_t filecount  = 1;                      // number of blocks

        // init the offsets etc
        const int block_offset_len = n_block_local + is_last;
        ullong*   block_offset     = reinterpret_cast<ullong*>(m_calloc(block_offset_len * sizeof(ullong)));
        const int block_center_len = n_block_local;
        int*      block_center     = reinterpret_cast<int*>(m_calloc(3 * block_center_len * sizeof(int)));
        const int block_box_len    = n_block_local;
        int*      block_box        = reinterpret_cast<int*>(m_calloc(6 * block_box_len * sizeof(int)));

        //................................................
        // do the data and get info for the offsets, centers and box
        for (iblock_t lbid = 0; lbid < n_block_local; ++lbid) {
            // get the associated GridBlock
            const iblock_t bid   = p4est_GetQuadIdOnLevel(mesh, il, lbid);
            p8est_tree_t*  tree  = p8est_tree_array_index(forest->trees, mesh->quad_to_tree[bid]);
            qdrt_t*        quad  = p8est_quadrant_array_index(&tree->quadrants, (bid - tree->quadrants_offset));
            GridBlock*     block = *(reinterpret_cast<GridBlock**>(quad->p.user_data));

            // move the file hyperslab forward
            const hsize_t fileoffset = (n_block_offset + lbid) * block_mem * field_lda;
            status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &fileoffset, &filestride, &filecount, &fileblock);
            m_assert(status >= 0, "failed to select the file hyperslab: block = %d, level = %d", lbid, il);

            // copy
            mem_ptr local_data = block->pointer(field) + local_mem_offset;
            status = H5Dwrite(dataset, M_HDF5_REAL, memdataspace, dataspace, dxpl, local_data);
            m_assert(status >= 0, "issue while writting the data: block = %d, level = %d", lbid, il);

            // store the offsets stuffs etc
            block_offset[lbid] = fileoffset;
            // it's node based so, the value is 1
            block_center[lbid * 3 + 0] = 1;
            block_center[lbid * 3 + 1] = 1;
            block_center[lbid * 3 + 2] = 1;
            // get the box of the block
            block_box[lbid * 6 + 0] = (int)(block->xyz(0) / block->hgrid(0)) - (int)((dump_ghost_) ? M_GS : 0);
            block_box[lbid * 6 + 1] = (int)(block->xyz(1) / block->hgrid(1)) - (int)((dump_ghost_) ? M_GS : 0);
            block_box[lbid * 6 + 2] = (int)(block->xyz(2) / block->hgrid(2)) - (int)((dump_ghost_) ? M_GS : 0);
            block_box[lbid * 6 + 3] = block_box[lbid * 6 + 0] + (int)(block_stride)-1;
            block_box[lbid * 6 + 4] = block_box[lbid * 6 + 1] + (int)(block_stride)-1;
            block_box[lbid * 6 + 5] = block_box[lbid * 6 + 2] + (int)(block_stride)-1;
        }
        if (is_last) {
            block_offset[n_block_local] = (n_block_offset + n_block_local) * block_mem * field_lda;
        }

        // create the datasets of the group
        hsize_t offset, local_size;
        hid_t   memspace;
        //................................................
        // the boxes: length = 1 babox_id per block
        string  bdsname("boxes");
        hsize_t box_len      = n_block_global;
        hid_t   boxdataspace = H5Screate_simple(1, &box_len, NULL);
        hid_t   boxdataset   = H5Dcreate(grp, bdsname.c_str(), babox_id, boxdataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        // create the mem dataspace and select the file hyperslab
        offset     = n_block_offset;
        local_size = block_box_len;
        memspace   = H5Screate_simple(1, &local_size, NULL);
        status     = H5Sselect_hyperslab(boxdataspace, H5S_SELECT_SET, &offset, NULL, &local_size, NULL);
        m_assert(status >= 0, "issue while selecting the hyperslab for the offset");
        status = H5Dwrite(boxdataset, babox_id, memspace, boxdataspace, dxpl, block_box);
        m_assert(status >= 0, "issue while writting offset");
        // close everything
        H5Sclose(memspace);
        H5Sclose(boxdataspace);
        H5Dclose(boxdataset);
        m_free(block_box);

        //................................................
        // the offsets - create the dataset/space
        string  odsname("data:offsets=0");
        hsize_t ofset_len       = n_block_global + 1;
        hid_t   offsetdataspace = H5Screate_simple(1, &ofset_len, NULL);
        hid_t   offsetdataset   = H5Dcreate(grp, odsname.c_str(), H5T_NATIVE_LLONG, offsetdataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        // create the mem dataspace and select the file hyperslab
        offset     = n_block_offset;
        local_size = block_offset_len;
        memspace   = H5Screate_simple(1, &local_size, NULL);
        status     = H5Sselect_hyperslab(offsetdataspace, H5S_SELECT_SET, &offset, NULL, &local_size, NULL);
        m_assert(status >= 0, "issue while selecting the hyperslab for the offset");
        status = H5Dwrite(offsetdataset, H5T_NATIVE_LLONG, memspace, offsetdataspace, dxpl, block_offset);
        m_assert(status >= 0, "issue while writting offset");
        // close everything
        H5Sclose(memspace);
        H5Sclose(offsetdataspace);
        H5Dclose(offsetdataset);
        m_free(block_offset);

        //................................................
        // the centers, defines if it's a cell (=0) or a node (=1) based info
        string  centername("boxcenter");
        hsize_t center_len      = n_block_global;
        hid_t   centerdataspace = H5Screate_simple(1, &center_len, NULL);
        hid_t   centerdataset   = H5Dcreate(grp, centername.c_str(), center_id, centerdataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        // create the mem dataspace and select the file hyperslab
        offset     = n_block_offset;
        local_size = block_center_len;
        memspace   = H5Screate_simple(1, &local_size, NULL);
        status     = H5Sselect_hyperslab(centerdataspace, H5S_SELECT_SET, &offset, NULL, &local_size, NULL);
        m_assert(status >= 0, "issue while selecting the hyperslab for the offset");
        status = H5Dwrite(centerdataset, center_id, memspace, centerdataspace, dxpl, block_center);
        m_assert(status >= 0, "issue while writting offset");
        // close everything
        H5Sclose(memspace);
        H5Sclose(centerdataspace);
        H5Dclose(centerdataset);
        m_free(block_center);

        //................................................
        // close everything for the current level
        H5Sclose(dataspace);
        H5Dclose(dataset);
        H5Gclose(grp);
    }

    //................................................
    // close the life
    H5Sclose(memdataspace);
    H5Tclose(center_id);
    H5Tclose(babox_id);
    H5Pclose(fapl);
    H5Pclose(dxpl);
    H5Fclose(hdf5_file);

    // reset the dump_ghost to false for next dump
    dump_ghost_ = false;
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief ask to dump the ghosts, this will be reset to false after the next dump
 * 
 * @param dump_ghost 
 */
void IOHDF5::dump_ghost(const bool dump_ghost) {
    dump_ghost_ = dump_ghost;
}

void IOHDF5::ApplyConstOpF(const qid_t* qid, GridBlock* block, const Field* fid) {
    //-------------------------------------------------------------------------
    hdf5_write_block_(qid, block, fid);
    //-------------------------------------------------------------------------
}

void IOHDF5::hdf5_write_header_(const ForestGrid* grid, const Field* field, const level_t min_lvl, const level_t max_lvl) {
    m_begin;
    //-------------------------------------------------------------------------
    // get some nice info
    rank_t rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    real_t L[3] = {grid->domain_length(0), grid->domain_length(1), grid->domain_length(2)};

    //................................................
    // only the master rank will create the file and write the metadata
    hid_t fapl, dxpl, fid, grp;
    if (rank == 0) {
        // Have only one rank to create and write metadata (header)
        fapl = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(fapl, MPI_COMM_SELF, MPI_INFO_NULL);
        H5Pset_coll_metadata_write(fapl, true);
        H5Pset_all_coll_metadata_ops(fapl, true);

        // Create the HDF5 file
        fid = H5Fcreate(filename_.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
        H5Pclose(fapl);
        m_assert(fid >= 0, "failed to create the hdf5 file");

        HDF5_WriteMetaData_(fid, field, min_lvl, max_lvl, L);
        // close the file
        H5Fclose(fid);
    }
    //-------------------------------------------------------------------------
    m_end;
}

void IOHDF5::HDF5_WriteMetaData_(hid_t fid, const Field* field, const level_t min_lvl, const level_t max_lvl, const real_t L[3]) {
    //-------------------------------------------------------------------------
    HDF5_AttributeString(fid, "version_name", "HyperClaw-V1.1");
    HDF5_AttributeString(fid, "plotfile_type", "VanillaHDF5");

    // setup the field name
    lid_t n_field = field->lda();
    HDF5_AttributeInt(fid, "num_components", 1, &n_field);
    char field_name[32], nickname[64];
    for (lda_t id = 0; id < n_field; ++id) {
        sprintf(field_name, "component_%d", id);
        sprintf(nickname, "%s_%d", field->name().c_str(), id);
        HDF5_AttributeString(fid, field_name, nickname);
    }

    // setup the dimensions
    lid_t ndim = 3;
    HDF5_AttributeInt(fid, "dim", 1, &ndim);

    // setup the current time
    real_t cur_time = 0.0;
    HDF5_AttributeReal(fid, "time", 1, &cur_time);

    // setup the level
    int finest = (int)(max_lvl - min_lvl);
    HDF5_AttributeInt(fid, "finest_level", 1, &finest);

    // setup the coordinate system, from AMReX_CoordSys 0 = cartesian, 2 = spherical
    int coord = 0;
    HDF5_AttributeInt(fid, "coordinate_system", 1, &coord);

    // setup the max level
    int nlevels = max_lvl - min_lvl + 1;
    HDF5_AttributeInt(fid, "num_levels", 1, &nlevels);

    // create the chombo local info
    hid_t grp = H5Gcreate(fid, "Chombo_global", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    HDF5_AttributeInt(grp, "SpaceDim", 1, &ndim);
    H5Gclose(grp);

    // get the bounding box type ready
    hid_t comp_dtype = H5Tcreate(H5T_COMPOUND, 6 * sizeof(int));
    H5Tinsert(comp_dtype, "lo_i", 0 * sizeof(int), H5T_NATIVE_INT);
    H5Tinsert(comp_dtype, "lo_j", 1 * sizeof(int), H5T_NATIVE_INT);
    H5Tinsert(comp_dtype, "lo_k", 2 * sizeof(int), H5T_NATIVE_INT);
    H5Tinsert(comp_dtype, "hi_i", 3 * sizeof(int), H5T_NATIVE_INT);
    H5Tinsert(comp_dtype, "hi_j", 4 * sizeof(int), H5T_NATIVE_INT);
    H5Tinsert(comp_dtype, "hi_k", 5 * sizeof(int), H5T_NATIVE_INT);

    // create the meta for the level groups
    for (level_t il = min_lvl; il <= max_lvl; ++il) {
        char level_name[32];
        sprintf(level_name, "level_%d", il - min_lvl);
        grp = H5Gcreate(fid, level_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        m_assert(grp >= 0, "error while writting in the hdf5 file");

        // the ref_ratio is given by the ratio between the # of point of 2 successive levels
        int ratio = m_min(max_lvl - il + 1, 2);  // (int)pow(2, max_lvl - il);
        HDF5_AttributeInt(grp, "ref_ratio", 1, &ratio);

        real_t cellsizes[3];
        for (lda_t id = 0; id < 3; ++id) {
            cellsizes[id] = 1.0 / ((double)((int)pow(2, il) * M_N));
        }
        // Visit has issues with vec_dx, and is ok with a single "dx" value
        HDF5_AttributeReal(grp, "Vec_dx", 3, cellsizes);
        HDF5_AttributeReal(grp, "dx", 1, &cellsizes[0]);  // For VisIt Chombo plot

        // get the box for the current level, at worse it's inside the bigbox
        real_t lo[3], hi[3];
        for (lda_t id = 0; id < 3; ++id) {
            lo[id] = 0.0;
            hi[id] = L[id];
        }
        HDF5_AttributeReal(grp, "prob_lo", 3, lo);
        HDF5_AttributeReal(grp, "prob_hi", 3, hi);
        // get the level resolution, at worse
        int domain[6];
        for (lda_t id = 0; id < 3; ++id) {
            domain[id]     = 0;
            domain[id + 3] = (int)(L[id] / cellsizes[id]) - 1;
        }

        hid_t aid         = H5Screate(H5S_SCALAR);
        hid_t domain_attr = H5Acreate(grp, "prob_domain", comp_dtype, aid, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(domain_attr, comp_dtype, domain);
        H5Aclose(domain_attr);
        H5Sclose(aid);

        // the type of the domain is cartesian  in every direction
        int type[3] = {0, 0, 0};
        HDF5_AttributeInt(grp, "domain_type", 3, type);

        // the level step is set to 0 (wtf is that?)
        int level_steps = 0;
        HDF5_AttributeInt(grp, "steps", 1, &level_steps);

        // again, wtf it that?
        int ngrow = 0;
        HDF5_AttributeInt(grp, "ngrow", 1, &ngrow);

        // we only IO one grid at the current time
        // int ngrid = 1;
        // HDF5_AttributeInt(grp, "ngrid", 1, &ngrid);
        HDF5_AttributeReal(grp, "time", 1, &cur_time);

        H5Gclose(grp);
    }

    H5Tclose(comp_dtype);

    //-------------------------------------------------------------------------
}

void IOHDF5::hdf5_write_block_(const qid_t* qid, GridBlock* block, const Field* fid) {
    //-------------------------------------------------------------------------
    //-------------------------------------------------------------------------
}