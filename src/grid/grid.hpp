#ifndef SRC_GRID_GRID_HPP_
#define SRC_GRID_GRID_HPP_

#include <p8est_mesh.h>

#include <iostream>
#include <limits>
#include <map>
#include <string>

#include "core/macros.hpp"
// #include "core/pointers.hpp"
#include "core/types.hpp"
#include "grid/field.hpp"
#include "grid/ghost.hpp"
#include "grid/gridblock.hpp"
#include "grid/gridcallback.hpp"
#include "operator/setvalues.hpp"
#include "tools/patch.hpp"
#include "tools/prof.hpp"
#include "wavelet/wavelet.hpp"

/**
 * @brief provide understandable grid management
 * 
 */
class Grid : public ForestGrid {
   protected:
    std::map<std::string, Field*> fields_;  //!< map of every field registered to this grid (the key is the field name, `field->name()`)

    Prof*    prof_   = nullptr;  //!< the profiler to use, may stay null
    Ghost*   ghost_  = nullptr;  //!< the ghost structure that handles one dimension of a field
    Wavelet* interp_ = nullptr;  //!< the Wavelet to use for all the multilevel operations

    bool   recursive_adapt_   = false;   //!< recursive adaptation or not
    real_t rtol_              = 1.0e-2;  //!< refinement tolerance, see @ref SetTol()
    real_t ctol_              = 1.0e-4;  //!< coarsening tolerance, see @ref SetTol()
    lid_t  n_quad_to_refine_  = 0;
    lid_t  n_quad_to_coarsen_ = 0;

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

    // void CopyFrom(const Grid* grid);
    void SetupMeshGhost();
    void DestroyMeshGhost();

    /**
     * @name Fields management
     * 
     * @{
     */
    [[nodiscard]] bidx_t NField() const { return static_cast<bidx_t>(fields_.size()); }

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
    // [[nodiscard]] inline lid_t NGhostFront() const {
    //     m_assert(interp_ != nullptr, "interp cannot be null");
    //     return interp_->nghost_front();
    // }
    // [[nodiscard]] inline lid_t NGhostBack() const {
    //     m_assert(interp_ != nullptr, "interp cannot be null");
    //     return interp_->nghost_back();
    // }
    inline void GhostLengthAdapt(bidx_t ghost_len[2]) const {
        ghost_len[0] = interp_->nghost_front();
        ghost_len[1] = interp_->nghost_back();
    }

    void GhostPull(Field* field, const BlockOperator* op) const;
    void GhostPull(Field* field, const bidx_t ghost_len_usr[2]) const;

    void GhostPull_SetLength(const Field* field, bidx_t ghost_len[2]) const;
    void GhostPull_Post(const Field* field, const sid_t ida, const bidx_t ghost_len[2]) const;
    void GhostPull_Wait(const Field* field, const sid_t ida, const bidx_t ghost_len[2]) const;
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

    [[nodiscard]] lid_t NQuadToAdapt() const { return (n_quad_to_refine_ + n_quad_to_coarsen_); }
    void                AddOneQuadToRefine() { n_quad_to_refine_ += 1; }
    void                AddOneQuadToCoarsen() { n_quad_to_coarsen_ += 8; }
    // void                AddOneQuadToAdapt() { ++n_quad_to_adapt_; }
    // void                AddQuadToAdapt(const sid_t n_quad) { n_quad_to_adapt_ += n_quad; }

    void SetTol(const real_t refine_tol, const real_t coarsen_tol);
    void SetRecursiveAdapt(const bool recursive_adapt) { recursive_adapt_ = recursive_adapt; }

    void StoreDetails(Field* criterion, Field* details);
    void MaxMinDetails(Field* criterion, real_t maxmin[2]);
    // void DistributionDetails(const iter_t id, const std::string folder, const std::string suffix, Field* criterion,
    //                          const short_t n_cat, const real_t max_category);
    // void DistributionDetailsInfiniteNorm(const iter_t id, const std::string folder, const std::string suffix, Field* criterion,
                            //  const short_t n_cat, const real_t max_category);

    void Refine(Field* field, const SetValue* expr = nullptr);
    void Coarsen(Field* field, const SetValue* expr = nullptr);
    void Adapt(Field* field, const SetValue* expr = nullptr);

    void Adapt(std::list<Patch>* patches);

    // void AdaptMagic(Field*  field, list<Patch*  > patches, cback_coarsen_citerion_t coarsen_crit, cback_refine_criterion_t refine_crit, void* criterion_ptr, cback_interpolate_t interp_fct, void* interp_ptr);
    void AdaptMagic(/* criterion */ Field* field_detail, std::list<Patch>* patches,
                    /* p4est coarsen/refine */ cback_coarsen_citerion_t coarsen_cback, cback_refine_criterion_t refine_cback, void* coarseref_cback_ptr,
                    /* p4est interpolate */ cback_interpolate_t interpolate_fct, void* interpolate_ptr);

    //    private:
    // void ExchangeStatus_PostStart_() const;
    // void ExchangeStatus_CompleteWait_() const;
    /**@}*/
};

#endif  // SRC_GRID_HPP_
