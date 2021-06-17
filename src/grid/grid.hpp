#ifndef SRC_GRID_GRID_HPP_
#define SRC_GRID_GRID_HPP_

#include <p8est_mesh.h>

#include <iostream>
#include <limits>
#include <map>
#include <string>

#include "core/macros.hpp"
#include "core/types.hpp"
#include "core/pointers.hpp"
\
#include "grid/gridblock.hpp"

#include "wavelet/wavelet.hpp"

#include "field.hpp"
#include "ghost.hpp"
#include "patch.hpp"
#include "prof.hpp"
#include "setvalues.hpp"
#include "gridcallback.hpp"

/**
 * @brief provide understandable grid management
 * 
 */
class Grid : public ForestGrid {
   protected:
    std::map<std::string, Field*  > fields_;  //!< map of every field registered to this grid (the key is the field name, `field->name()`)

    Prof*    prof_   = nullptr;  //!< the profiler to use, may stay null
    Ghost*   ghost_  = nullptr;  //!< the ghost structure that handles one dimension of a field
    Wavelet* interp_ = nullptr;  //!< the Wavelet to use for all the multilevel operations

    bool   recursive_adapt_ = false;   //!< recursive adaptation or not
    real_t rtol_            = 1.0e-2;  //!< refinement tolerance, see @ref SetTol()
    real_t ctol_            = 1.0e-4;  //!< coarsening tolerance, see @ref SetTol()
    lid_t  n_quad_to_adapt_ = 0;

    level_t level_limit_max_ = P8EST_QMAXLEVEL;  //!< max level of a quadrant in the mesh
    level_t level_limit_min_ = 0;                //!< min level of a quadrant

    void* cback_criterion_ptr_   = nullptr;  //!< temporary pointer to be used in the criterion callback functions
    void* cback_interpolate_ptr_ = nullptr;  //!< temporary pointer to be used in the interpolation callback functions

    // MPI_Win neighbor_status_window_ = MPI_WIN_NULL;  //!< window to access the status of ghost blocks
    // short*  neighbor_status_        = nullptr;       //!< status of every block: true = coarsen, false = do not coarsen

   public:
    explicit Grid();
    Grid(const level_t ilvl, const bool isper[3], const lid_t l[3], MPI_Comm comm, Prof* const prof);
    ~Grid();

    [[nodiscard]] size_t LocalMemSize() const;
    [[nodiscard]] size_t LocalNumDof() const;
    [[nodiscard]] size_t GlobalNumDof() const;

    [[nodiscard]] level_t level_limit_max() const { return level_limit_max_; }
    [[nodiscard]] level_t level_limit_min() const { return level_limit_min_; }
    void                  level_limit(const level_t min, const level_t max);

    [[nodiscard]] Wavelet* interp() const { return interp_; }
    [[nodiscard]] Prof*    profiler() const { return prof_; }

    bool HasProfiler() { return prof_ == nullptr; }

    void CopyFrom(const Grid* grid);
    void SetupMeshGhost();
    void DestroyMeshGhost();

    /**
     * @name Fields management
     * 
     * @{
     */
    [[nodiscard]] lid_t NField() const { return static_cast<lid_t>(fields_.size()); }

    [[nodiscard]] auto FieldBegin() const { return fields_.cbegin(); }
    [[nodiscard]] auto FieldEnd() const { return fields_.cend(); }

    bool IsAField(const Field* field) const;
    void AddField(Field* field);
    void DeleteField(const Field* field);
    void ResetFields(const std::map<std::string, Field*>* fields);

    /**@}*/

    /**
     * @name Ghost management
     * 
     * @{
     */
    [[nodiscard]] inline lid_t NGhostFront() const {
        m_assert(interp_ != nullptr, "interp cannot be null");
        return interp_->nghost_front();
    }
    [[nodiscard]] inline lid_t NGhostBack() const {
        m_assert(interp_ != nullptr, "interp cannot be null");
        return interp_->nghost_back();
    }
    void GhostPull(Field* field) const;
    void GhostPull_Post(const Field* field, const sid_t ida) const;
    void GhostPull_Wait(const Field* field, const sid_t ida) const;
    /**@}*/

    /**
     * @name Grid adaptation
     * 
     * @{
     */
    [[nodiscard]] real_t rtol() const { return rtol_; }
    [[nodiscard]] real_t ctol() const { return ctol_; }
    [[nodiscard]] bool   recursive_adapt() const { return recursive_adapt_; }
    [[nodiscard]] void*  cback_criterion_ptr() const { return cback_criterion_ptr_; }
    [[nodiscard]] void*  cback_interpolate_ptr() const { return cback_interpolate_ptr_; }

    [[nodiscard]] lid_t n_quad_to_adapt() const { return n_quad_to_adapt_; }
    void  AddOneQuadToAdapt() { ++n_quad_to_adapt_; }
    void  AddQuadToAdapt(const sid_t n_quad) { n_quad_to_adapt_ += n_quad; }

    void SetTol(const real_t refine_tol, const real_t coarsen_tol);
    void SetRecursiveAdapt(const bool recursive_adapt) { recursive_adapt_ = recursive_adapt; }

    void Refine(Field*  field);
    void Coarsen(Field*  field);
    void StoreDetails(Field*  criterion, Field*  details);

    void Adapt(Field*  field);
    void Adapt(Field*  field, SetValue*  expression);
    void Adapt(std::list<Patch>* patches);

    // void AdaptMagic(Field*  field, list<Patch*  > patches, cback_coarsen_citerion_t coarsen_crit, cback_refine_criterion_t refine_crit, void* criterion_ptr, cback_interpolate_t interp_fct, void* interp_ptr);
    void AdaptMagic(/* criterion */ Field*  field_detail, std::list<Patch>* patches,
                    /* p4est coarsen/refine */ cback_coarsen_citerion_t coarsen_cback, cback_refine_criterion_t refine_cback, void* coarseref_cback_ptr,
                    /* p4est interpolate */ cback_interpolate_t interpolate_fct, void* interpolate_ptr);

//    private:
    // void ExchangeStatus_PostStart_() const;
    // void ExchangeStatus_CompleteWait_() const;
    /**@}*/
};

#endif  // SRC_GRID_HPP_
