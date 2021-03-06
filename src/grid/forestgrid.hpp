#ifndef SRC_FORESTGRID_HPP_
#define SRC_FORESTGRID_HPP_

#include <string>

#include "core/macros.hpp"
#include "core/types.hpp"
#include "tools/toolsp4est.hpp"

/**
 * @brief handles the p4est structures. This class is a minimal grid, intended to be used by every other class
 * that needs to access the p4est structures (@ref Ghost, @ref Partitioner, ...)
 * 
 */
class ForestGrid {
   protected:
    level_t global_min_level_ = -1;      //!< store the minimum level on the grid
    level_t global_max_level_ = -1;      //!< store the minimum level on the grid
    real_t  domain_length_[3];    //!< length of the domain, given by the tree aspect ratio
    bool    domain_periodic_[3];  //!< boolean indicating if a direction is periodic

    p8est_t*       p4est_forest_ = nullptr;  //!< the p8est_t structure
    p8est_mesh_t*  p4est_mesh_   = nullptr;  //!< the p8est_mesh_t structure
    p8est_ghost_t* p4est_ghost_  = nullptr;  //!< the p8est_ghost_t structure

    bool is_connect_owned_ = false; /**<indicate if we have to free the connectivity */
    bool is_mesh_valid_    = false; /**<indicate that the mesh and the ghost structures are up-to-date */

    BlockDataType block_type_ = M_NULLTYPE;

   public:
    /**
     * @name constructors and destructor
     * 
     * @{
     */
    explicit ForestGrid();
    ForestGrid(const level_t ilvl, const bool isper[3], const lid_t l[3], const BlockDataType block_type, MPI_Comm comm);
    ~ForestGrid();
    /** @} */

    // void CopyFrom(const ForestGrid* grid);

    [[nodiscard]] BlockDataType    block_type() const { return block_type_; }
    [[nodiscard]] real_t           domain_length(const sid_t id) const { return domain_length_[id]; }
    [[nodiscard]] bool             domain_periodic(const sid_t id) const { return domain_periodic_[id]; }
    [[nodiscard]] inline bool      is_mesh_valid() const { return is_mesh_valid_; }
    [[nodiscard]] p4est_Essentials p4estEssentials() const;

    /**
     * @name p4est managmeent
     * 
     * @{
     */
    [[nodiscard]] inline int                   mpirank() const { return p4est_forest_->mpirank; }
    [[nodiscard]] inline int                   mpisize() const { return p4est_forest_->mpisize; }
    [[nodiscard]] inline MPI_Comm              mpicomm() const { return p4est_forest_->mpicomm; }
    [[nodiscard]] inline p8est_t*              p4est_forest() const { return p4est_forest_; }
    [[nodiscard]] inline p8est_mesh_t*         p4est_mesh() const { return p4est_mesh_; }
    [[nodiscard]] inline p8est_ghost_t*        p4est_ghost() const { return p4est_ghost_; }
    [[nodiscard]] inline p8est_connectivity_t* p4est_connect() const { return p4est_forest_->connectivity; }
    [[nodiscard]] inline bidx_t                local_num_quadrants() const { return p4est_forest_->local_num_quadrants; }
    [[nodiscard]] inline long                  global_num_quadrants() const { return p4est_forest_->global_num_quadrants; }

    [[nodiscard]] level_t MaxLevel() const;
    [[nodiscard]] level_t MinLevel() const;
    [[nodiscard]] real_t  FinestH() const;
    [[nodiscard]] real_t  CoarsestH() const;

    void DumpLevels(const iter_t id, const std::string folder, const std::string suffix = "") const;
    /** @} */

    /**
     * @name Mesh and Ghost structure management
     * @{
     */
    void DestroyP4estMeshAndGhost();
    void SetupP4estMeshAndGhost();
    /** @} */
};

#endif  // SRC_FORESTGRID_HPP_
