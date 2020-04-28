#ifndef SRC_FORESTGRID_HPP_
#define SRC_FORESTGRID_HPP_

#include <p8est.h>
#include <p8est_extended.h>

#include "murphy.hpp"

/**
 * @brief handles the p4est structures. This class is a minimal grid, intended to be used by every other class
 * that needs to access the p4est structures (@ref Ghost, @ref Partitioner, ...)
 * 
 */
class ForestGrid {
   protected:
    real_t domain_length_[3];    //!< length of the domain, given by the tree aspect ratio
    bool   domain_periodic_[3];  //!< boolean indicating if a direction is periodic

    p8est_t*       forest_ = nullptr;  //!< the p8est_t structure
    p8est_mesh_t*  mesh_   = nullptr;  //!< the p8est_mesh_t structure
    p8est_ghost_t* ghost_  = nullptr;  //!< the p8est_ghost_t structure

    bool is_mesh_valid_ = false; /**<indicate that the mesh and the ghost structures are up-to-date */

   public:
    ForestGrid(const lid_t ilvl, const bool isper[3], const lid_t l[3], const size_t datasize, MPI_Comm comm);
    ~ForestGrid();

    real_t      domain_length(const sid_t id) { return domain_length_[id]; }
    bool        domain_periodic(const sid_t id) { return domain_periodic_[id]; }
    inline bool is_mesh_valid() const { return is_mesh_valid_; }

    /**
    * @name p4est managmeent
    * 
    * @{
    */
    inline int                   mpirank() const { return forest_->mpirank; }
    inline int                   mpisize() const { return forest_->mpisize; }
    inline MPI_Comm              mpicomm() const { return forest_->mpicomm; }
    inline p8est_t*              forest() const { return forest_; }
    inline p8est_mesh_t*         mesh() const { return mesh_; }
    inline p8est_ghost_t*        ghost() const { return ghost_; }
    inline p8est_connectivity_t* connect() const { return forest_->connectivity; }

    inline lid_t local_num_quadrants() const { return forest_->local_num_quadrants; }
    /** @} */

    /**
     * @name Mesh and Ghost structure management
     * @{
     */
    void ResetP4estGhostMesh();
    void SetupP4estGhostMesh();
    /** @} */
};

#endif  // SRC_FORESTGRID_HPP_
