#ifndef SRC_FOREST_GRID_HPP_
#define SRC_FOREST_GRID_HPP_

#include <p8est.h>
#include <p8est_extended.h>

#include "murphy.hpp"

class ForestGrid {
   protected:
    real_t domain_length_[3];
    bool   domain_periodic_[3];

    p8est_t*       forest_ = nullptr;
    p8est_mesh_t*  mesh_   = nullptr;
    p8est_ghost_t* ghost_  = nullptr;

    bool is_mesh_valid_ = false; /**<indicate that the mesh it no usable */

   public:
    inline int                   mpirank() const { return forest_->mpirank; }
    inline int                   mpisize() const { return forest_->mpisize; }
    inline MPI_Comm              mpicomm() const { return forest_->mpicomm; }
    inline p8est_t*              forest() const { return forest_; }
    inline p8est_mesh_t*         mesh() const { return mesh_; }
    inline p8est_ghost_t*        ghost() const { return ghost_; }
    inline p8est_connectivity_t* connect() const { return forest_->connectivity; }

    real_t domain_length(const sid_t id) { return domain_length_[id]; }
    bool   domain_periodic(const sid_t id) { return domain_periodic_[id]; }

    inline bool  is_mesh_valid() const { return is_mesh_valid_; }
    inline lid_t local_num_quadrants() const { return forest_->local_num_quadrants; }

    ForestGrid(const lid_t ilvl, const bool isper[3], const lid_t l[3], const size_t datasize, MPI_Comm comm);
    ~ForestGrid();

    /**
     * @name Mesh and Ghost management
     * @{
     */
    void ResetP4estGhostMesh();
    void SetupP4estGhostMesh();
    /** @} */
};

#endif  // SRC_FOREST_GRID_HPP_
