#ifndef SRC_GRID_HPP_
#define SRC_GRID_HPP_

#include <p8est_mesh.h>

#include <iostream>
#include <map>
#include <string>

#include "field.hpp"
#include "murphy.hpp"
#include "prof.hpp"

/**
 * @brief Multiresolution grid structure that wraps the forest from p4est, the mesh and the ghost object
 * 
 * The 3 objects are linked together since we need a forest and a ghost to have a mesh
 * 
 */
class Grid {
   protected:
    bool is_mesh_valid_ = false; /**<indicate that the mesh it no usable */

    p8est_t*       forest_ = nullptr; /**< the p8est_t topology */
    p8est_ghost_t* ghost_  = nullptr; /**< the p8est_ghost data-structure */
    p8est_mesh_t*  mesh_   = nullptr; /**< the p8est_mesh access */

    ldamap_t lda_; /*<store the dimension of every field, using its tag */

    Prof* prof_ = nullptr; /**< the profiler */

   public:
    Grid(const lid_t ilvl, const bool isper[3], const lid_t l[3], MPI_Comm comm, Prof* prof);
    ~Grid();

    p8est_t*              forest() { return forest_; }
    p8est_mesh_t*         mesh() { return mesh_; }
    p8est_ghost_t*        ghost() { return ghost_; }
    p8est_connectivity_t* connect() { return forest_->connectivity; }

    bool     is_mesh_valid() const { return is_mesh_valid_; }
    int      mpirank() const { return forest_->mpirank; }
    MPI_Comm mpicomm() const { return forest_->mpicomm; }
    lid_t    local_num_quadrants() const { return forest_->local_num_quadrants; }

    size_t LocalMemSize() const;
    size_t LocalNumDof() const;
    size_t GlobalNumDof() const;
    bool   IsAField(const Field* field) const;

    void AddField(Field* field);
    void DeleteField(Field* field);
};

#endif  // SRC_GRID_HPP_