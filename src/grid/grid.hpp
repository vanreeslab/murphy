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

#include "grid/gridblock.hpp"

#include "wavelet/wavelet.hpp"

#include "field.hpp"
#include "ghost.hpp"
#include "patch.hpp"
#include "prof.hpp"
#include "setvalues.hpp"
#include "gridcallback.hpp"

/**
 * @brief implements the grid management and the related responsabilities on top of ForestGrid
 * 
 */
class Grid : public ForestGrid {
   protected:
    std::map<std::string, m_ptr<Field> > fields_;  //!< map of every field registered to this grid (the key is the field name, `field->name()`)

    Prof*    prof_   = nullptr;  //!< the profiler to use, may stay null
    Ghost*   ghost_  = nullptr;  //!< the ghost structure that handles one dimension of a field
    Wavelet* interp_ = nullptr;  //!< the Wavelet to use for all the multilevel operations

    bool   recursive_adapt_ = false;   //!< recursive adaptation or not
    real_t rtol_            = 1.0e-2;  //!< refinement tolerance, see @ref SetTol()
    real_t ctol_            = 1.0e-4;  //!< coarsening tolerance, see @ref SetTol()
    lid_t  n_quad_to_adapt_ = 0;

    void* cback_criterion_ptr_   = nullptr;  //!< temporary pointer to be used in the criterion callback functions
    void* cback_interpolate_ptr_ = nullptr;  //!< temporary pointer to be used in the interpolation callback functions

   public:
    explicit Grid();
    Grid(const lid_t ilvl, const bool isper[3], const lid_t l[3], MPI_Comm comm, Prof* prof);
    ~Grid();

    size_t LocalMemSize() const;
    size_t LocalNumDof() const;
    size_t GlobalNumDof() const;

    Wavelet* interp() const { return interp_; }
    Prof*    profiler() const { return prof_; }

    bool HasProfiler() { return prof_ == nullptr; }

    void CopyFrom(m_ptr<const Grid> grid);
    void SetupGhost();
    void DestroyGhost();

    /**
     * @name Fields management
     * 
     * @{
     */
    lid_t NField() const { return (lid_t)(fields_.size()); }

    auto FieldBegin() const { return fields_.cbegin(); }
    auto FieldEnd() const { return fields_.cend(); }

    bool IsAField(m_ptr<const Field> field) const;
    void AddField(m_ptr<Field> field);
    void DeleteField(m_ptr<const Field> field);
    void ResetFields(m_ptr<const std::map<std::string, m_ptr<Field> > > fields);

    /**@}*/

    /**
     * @name Ghost management
     * 
     * @{
     */
    inline lid_t NGhostFront() const {
        m_assert(interp_ == nullptr, "interp cannot be null");
        return interp_->nghost_front();
    }
    inline lid_t NGhostBack() const {
        m_assert(interp_ == nullptr, "interp cannot be null");
        return interp_->nghost_back();
    }
    void GhostPull(m_ptr<Field> field);
    void GhostPull_Post(m_ptr<const Field> field, const sid_t ida);
    void GhostPull_Wait(m_ptr<const Field> field, const sid_t ida);
    /**@}*/

    /**
     * @name Grid adaptation
     * 
     * @{
     */
    real_t rtol() const { return rtol_; }
    real_t ctol() const { return ctol_; }
    bool   recursive_adapt() const { return recursive_adapt_; }
    void*  cback_criterion_ptr() const { return cback_criterion_ptr_; }
    void*  cback_interpolate_ptr() const { return cback_interpolate_ptr_; }

    lid_t n_quad_to_adapt() const { return n_quad_to_adapt_; }
    void  AddOneQuadToAdapt() { ++n_quad_to_adapt_; }
    void  AddQuadToAdapt(const sid_t n_quad) { n_quad_to_adapt_ += n_quad; }

    void SetTol(const real_t refine_tol, const real_t coarsen_tol);
    void SetRecursiveAdapt(const bool recursive_adapt) { recursive_adapt_ = recursive_adapt; }

    void Refine(m_ptr<Field> field);
    void Coarsen(m_ptr<Field> field);
    void Adapt(m_ptr<Field> field);
    void Adapt(m_ptr<Field> field, m_ptr<SetValue> expression);
    void Adapt(m_ptr<std::list<Patch> > patches);
    void Adapt(m_ptr<const Field> field, cback_coarsen_citerion_t coarsen_crit, cback_refine_criterion_t refine_crit, void* criterion_ptr, cback_interpolate_t interp, void* interp_ptr);
    void DumpDetails(m_ptr<Field> criterion, m_ptr<Field> details);
    /**@}*/
};

#endif  // SRC_GRID_HPP_
