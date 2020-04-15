#include "forestgrid.hpp"

ForestGrid::ForestGrid(const lid_t ilvl, const bool isper[3], const lid_t l[3], const size_t datasize, MPI_Comm comm) {
    m_begin;
    m_assert(ilvl >= 0, "the init level has to be >= 0");
    m_assert(ilvl <= P8EST_MAXLEVEL, "the init level has to be <= P8EST_MAXLEVEL");
    //-------------------------------------------------------------------------
    // create the connect as a box of L[0]xL[1]xL[2] trees
    p8est_connectivity_t* connect = p8est_connectivity_new_brick(l[0], l[1], l[2], isper[0], isper[1], isper[2]);
    // create the forest at a given level, the associated ghost and mesh object
    forest_ = p8est_new_ext(comm, connect, 0, ilvl, 1, datasize, nullptr, nullptr);
    // set the pointer to null
    forest_->user_pointer = nullptr;
    // set the ghost and the mesh
    SetupP4estGhostMesh();
    //-------------------------------------------------------------------------
    m_end;
}


ForestGrid::~ForestGrid(){
    m_begin;
    //-------------------------------------------------------------------------
    // delete the ghost and the mesh
    ResetP4estGhostMesh();
    // destroy the connectivity and the forest
    p8est_connectivity_t* connect = forest_->connectivity;
    p8est_destroy(forest_);
    p8est_connectivity_destroy(connect);
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief destroys the mesh and the ghost p4est structure and unvalidate the mesh
 * 
 */
void ForestGrid::ResetP4estGhostMesh() {
    m_begin;
    //-------------------------------------------------------------------------
    // destroy the structures
    p8est_mesh_destroy(mesh_);
    p8est_ghost_destroy(ghost_);
    // unvalidate the mesh
    is_mesh_valid_ = false;
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief reinitializes the mesh and the ghost p4est structure and validate the mesh
 * 
 */
void ForestGrid::SetupP4estGhostMesh() {
    m_begin;
    //-------------------------------------------------------------------------
    ghost_ = p8est_ghost_new(forest_, P8EST_CONNECT_FULL);
    mesh_  = p8est_mesh_new_ext(forest_, ghost_, 1, 1, P8EST_CONNECT_FULL);
    
    is_mesh_valid_ = true;
    //-------------------------------------------------------------------------
    m_end;
}