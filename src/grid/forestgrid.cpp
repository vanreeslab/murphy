#include "forestgrid.hpp"

#include <p8est_extended.h>
#include "tools/toolsblocktypes.hpp"

/**
 * @brief Construct an empty ForestGrid
 * 
 */
ForestGrid::ForestGrid() {
    m_begin;
    //-------------------------------------------------------------------------
    p4est_forest_     = nullptr;
    p4est_mesh_       = nullptr;
    p4est_ghost_      = nullptr;
    is_mesh_valid_    = false;
    is_connect_owned_ = false;
    block_type_       = M_NULLTYPE;
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief Construct a new Forest Grid, initialize the p4est structures.
 * The forest is a uniform resolution forest, composed of l[0]xl[1]xl[2] octrees at refinement level ilvl
 * 
 * @param ilvl the constant refinement level in one tree
 * @param isper the periodicity of the whole domain, per direction
 * @param l the aspect ratio of the domain, in trees
 * @param datasize the datasize to give to p4est
 * @param comm the communicator to build the forest
 */
ForestGrid::ForestGrid(const level_t ilvl, const bool isper[3], const lid_t l[3], const BlockDataType block_type, MPI_Comm comm) {
    m_begin;
    m_assert(ilvl >= 0, "the init level has to be >= 0");
    m_assert(ilvl <= P8EST_MAXLEVEL, "the init level has to be <= P8EST_MAXLEVEL");
    //-------------------------------------------------------------------------
    // store the type of block managed by the grid
    block_type_ = block_type;
    // create the connect as a box of L[0]xL[1]xL[2] trees
    is_connect_owned_             = true;
    p8est_connectivity_t* connect = p8est_connectivity_new_brick(l[0], l[1], l[2], isper[0], isper[1], isper[2]);
    // create the forest at a given level, the associated ghost and mesh object
    int comm_size;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    // int min_quad  = ((int)pow(2, 3 * ilvl)) / comm_size;
    // p4est_forest_ = p8est_new_ext(comm, connect, min_quad, 0, 1, datasize, nullptr, nullptr);
    p4est_forest_ = p8est_new_ext(comm, connect,0, ilvl, 1, BlockPointerSize(block_type), nullptr, nullptr);
    // forest_ - p8est_new(comm,connect,datasize,nullptr,nullptr);
    // set the pointer to null
    p4est_forest_->user_pointer = nullptr;
    // store the domain periodicity and domain size
    domain_length_[0]   = (real_t)l[0];
    domain_length_[1]   = (real_t)l[1];
    domain_length_[2]   = (real_t)l[2];
    domain_periodic_[0] = isper[0];
    domain_periodic_[1] = isper[1];
    domain_periodic_[2] = isper[2];
    //-------------------------------------------------------------------------
    m_assert(p4est_forest_->global_num_quadrants == ((int)pow(2, 3 * ilvl) * l[0] * l[1] * l[2]), "the number of global quad = %ld must match the initial level number = %d", p4est_forest_->global_num_quadrants, (int)pow(2, 3 * ilvl));
    m_end;
}

// /**
//  * @brief copies the p4est forest of another grid while keeping the block address
//  * 
//  * @param grid the other grid to copy
//  */
// void ForestGrid::CopyFrom(const ForestGrid*  grid) {
//     m_begin;
//     //-------------------------------------------------------------------------
//     // copy the existing forest, including the memory adress to the GridBlock
//     is_connect_owned_ = false;
//     p4est_forest_     = p8est_copy(grid->p4est_forest(), 1);
//     // set the pointer to null
//     p4est_forest_->user_pointer = nullptr;
//     // store the domain periodicity and domain size
//     domain_length_[0]   = grid->domain_length(0);
//     domain_length_[1]   = grid->domain_length(1);
//     domain_length_[2]   = grid->domain_length(2);
//     domain_periodic_[0] = grid->domain_periodic(0);
//     domain_periodic_[1] = grid->domain_periodic(1);
//     domain_periodic_[2] = grid->domain_periodic(2);
//     //-------------------------------------------------------------------------
//     m_end;
// }

/**
 * @brief Destroy the Forest Grid, reset the ghost and the mesh and clean everything from the p4est side
 * 
 */
ForestGrid::~ForestGrid() {
    m_begin;
    //-------------------------------------------------------------------------
    // delete the ghost and the mesh
    DestroyP4estMeshAndGhost();
    // destroy the connectivity and the forest
    p8est_connectivity_t* connect = nullptr;
    if (is_connect_owned_) {
        connect = p4est_forest_->connectivity;
    }
    p8est_destroy(p4est_forest_);
    if (is_connect_owned_) {
        p8est_connectivity_destroy(connect);
    }
    //--------------------------------------------------------------------------
    m_end;
}

/**
 * @brief returns the essential information about p4est, cfr p4est_Essentials
 * 
 * @return p4est_Essentials 
 */
p4est_Essentials ForestGrid::p4estEssentials() const {
    m_begin;
    //--------------------------------------------------------------------------
    p4est_Essentials essential{domain_periodic_, p4est_forest_, p4est_mesh_, p4est_ghost_};
    return essential;
    //--------------------------------------------------------------------------
    m_end;
}

/**
 * @brief destroys the mesh and the ghost p4est structure and unvalidate the mesh.
 * 
 * This is a mandatory step when changing the grid
 * 
 */
void ForestGrid::DestroyP4estMeshAndGhost() {
    m_begin;
    //-------------------------------------------------------------------------
    // destroy the structures
    if (p4est_mesh_ != nullptr) {
        p8est_mesh_destroy(p4est_mesh_);
        p4est_mesh_ = nullptr;
        m_verb("destroyed the old mesh!!");
    }
    if (p4est_ghost_ != nullptr) {
        p8est_ghost_destroy(p4est_ghost_);
        p4est_ghost_ = nullptr;
        m_verb("destroyed the old ghosts!!");
    }
    // invalidate the levels
    global_max_level_ = -1;
    global_min_level_ = -1;

    // unvalidate the mesh
    is_mesh_valid_ = false;
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief reinitializes the mesh and the ghost p4est structure and validate the mesh
 * 
 * This a mandatory step when changing the grid
 * 
 */
void ForestGrid::SetupP4estMeshAndGhost() {
    m_begin;
    m_assert(p4est_ghost_ == nullptr && p4est_mesh_ == nullptr, "cannot initialize something that already exist");
    //-------------------------------------------------------------------------
    p4est_ghost_ = p8est_ghost_new(p4est_forest_, P8EST_CONNECT_FULL);
    m_verb("creating a new mesh");
    p4est_mesh_ = p8est_mesh_new_ext(p4est_forest_, p4est_ghost_, 1, 1, P8EST_CONNECT_FULL);
    m_verb("done with mesh");
    is_mesh_valid_ = true;

    {  // compute the levels: MAX
        level_t il = P8EST_QMAXLEVEL;
        for (il = P8EST_QMAXLEVEL; il >= 0; --il) {
            if (p4est_NumQuadOnLevel(p4est_mesh_, il) != 0) {
                break;
            }
        }
        m_assert(sizeof(level_t) == sizeof(short), "the MPI call is done for a char");
        MPI_Allreduce(&il, &global_max_level_, 1, MPI_SHORT, MPI_MAX, MPI_COMM_WORLD);
        m_assert(global_max_level_ >= 0, "the level=%d must be >=0", global_max_level_);
    }

    {  // compute the levels: MIN
        level_t il = 0;
        for (il = 0; il <= P8EST_QMAXLEVEL; ++il) {
            if (p4est_NumQuadOnLevel(p4est_mesh_, il) != 0) {
                break;
            }
        }

        // get the max among the cpus
        m_assert(sizeof(level_t) == sizeof(short), "the MPI call is done for a char");
        MPI_Allreduce(&il, &global_min_level_, 1, MPI_SHORT, MPI_MIN, MPI_COMM_WORLD);
        m_assert(global_min_level_ >= 0, "the level=%d must be >=0", global_min_level_);
    }

    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief returns the max active level on the current grid
 * 
 * @return level_t 
 */
level_t ForestGrid::MaxLevel() const {
    m_assert(is_mesh_valid(), "the mesh must be valid to return the max level");
    m_assert(global_max_level_ > -1, "the max level has not been computed yet");
    return global_max_level_;
    // //-------------------------------------------------------------------------
    // level_t il = P8EST_QMAXLEVEL;
    // for (il = P8EST_QMAXLEVEL; il >= 0; --il) {
    //     if (p4est_NumQuadOnLevel(p4est_mesh_, il) != 0) {
    //         break;
    //     }
    // }
    // level_t global_level = P8EST_QMAXLEVEL;
    // m_assert(sizeof(level_t) == sizeof(short), "the MPI call is done for a char");
    // MPI_Allreduce(&il, &global_level, 1, MPI_SHORT, MPI_MAX, MPI_COMM_WORLD);
    // m_assert(global_level >= 0, "the level=%d must be >=0", global_level);
    // return global_level;
    // //-------------------------------------------------------------------------
}

/**
 * @brief returns the min active level on the current grid
 * 
 * @return level_t 
 */
level_t ForestGrid::MinLevel() const {
    m_assert(is_mesh_valid(), "the mesh must be valid to return the max level");
    m_assert(global_min_level_ > -1, "the max level has not been computed yet");
    return global_min_level_;
    // //-------------------------------------------------------------------------
    // level_t il = 0;
    // for (il = 0; il <= P8EST_QMAXLEVEL; ++il) {
    //     if (p4est_NumQuadOnLevel(p4est_mesh_, il) != 0) {
    //         break;
    //     }
    // }

    // // get the max among the cpus
    // level_t global_level = 0;
    // m_assert(sizeof(level_t) == sizeof(short),"the MPI call is done for a char");
    // MPI_Allreduce(&il, &global_level, 1, MPI_SHORT, MPI_MIN, MPI_COMM_WORLD);
    // m_assert(global_level >= 0, "the level=%d must be >=0", global_level);
    // return global_level;
    // //-------------------------------------------------------------------------
}

real_t ForestGrid::FinestH() const {
    // m_begin;
    //-------------------------------------------------------------------------
    level_t max_level = MaxLevel();
    m_assert(max_level >= 0, "the level=%d must be >=0", max_level);
    return p4est_QuadLen(max_level) / ((real_t)M_N);
    //-------------------------------------------------------------------------
    // m_end;
}
real_t ForestGrid::CoarsestH() const {
    // m_begin;
    //-------------------------------------------------------------------------
    level_t min_level = MinLevel();
    m_assert(min_level >= 0, "the level=%d must be >=0", min_level);
    return p4est_QuadLen(min_level)/((real_t)M_N);
    //-------------------------------------------------------------------------
    // m_end;
}

/**
 * @brief Dump the number of blocks on each level
 * 
 * @param id the id put in the begining of the IO line
 */
void ForestGrid::DumpLevels(const iter_t id, const std::string folder, const std::string suffix) const {
    m_assert(is_mesh_valid(), "the mesh must be valid to return the max level");
    //-------------------------------------------------------------------------
    long_t n_per_level[P8EST_MAXLEVEL];
    long_t local_n_per_level[P8EST_MAXLEVEL];
    for (level_t il = 0; il < P8EST_MAXLEVEL; ++il) {
        local_n_per_level[il] = p4est_NumQuadOnLevel(p4est_mesh_, il);
    }
    // get the sum over the comm world
    MPI_Reduce(local_n_per_level, n_per_level, P8EST_MAXLEVEL, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    // get the rank
    rank_t rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // dump it
    FILE* file_diag;
    if (rank == 0) {
        file_diag = fopen(std::string(folder + "/grid_levels" + suffix + ".data").c_str(), "a+");
        fprintf(file_diag, "%d", id);
        for (level_t il = 0; il < P8EST_MAXLEVEL; ++il) {
            fprintf(file_diag, ";%ld", n_per_level[il]);
        }
        fprintf(file_diag, "\n");
        fclose(file_diag);
    }

    //-------------------------------------------------------------------------
}
