#ifndef SRC_GRID_HPP_
#define SRC_GRID_HPP_

#include <p8est_mesh.h>

#include <iostream>
#include <limits>
#include <unordered_map>
#include <string>

#include "field.hpp"
#include "ghost.hpp"
#include "gridblock.hpp"
#include "interpolator.hpp"
#include "murphy.hpp"
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
    std::unordered_map<std::string, Field*> fields_;  //!< map of every field registered to this grid (the key is the field name, `field->name()`)

    Prof*                 prof_   = nullptr;  //!< the profiler to use, may stay null
    Ghost*                ghost_  = nullptr;  //!< the ghost structure that handles one dimension of a field
    InterpolatingWavelet* interp_ = nullptr;  //!< the interpolator to use for all the multilevel operations

    bool   recursive_adapt_ = false;   //!< recursive adaptation or not
    real_t rtol_            = 1.0e-2;  //!< refinement tolerance, see @ref SetTol()
    real_t ctol_            = 1.0e-4;  //!< coarsening tolerance, see @ref SetTol()

    void* cback_criterion_ptr_ = nullptr;  //!< temporary pointer to be used in the criterion callback functions
    void* cback_interpolate_ptr_ = nullptr;  //!< temporary pointer to be used in the interpolation callback functions

   public:
    explicit Grid();
    Grid(const lid_t ilvl, const bool isper[3], const lid_t l[3], MPI_Comm comm, Prof* prof);
    ~Grid();

    size_t LocalMemSize() const;
    size_t LocalNumDof() const;
    size_t GlobalNumDof() const;

    InterpolatingWavelet* interp() { return interp_; }

    Prof* profiler() { return prof_; }
    bool  HasProfiler() { return prof_ != nullptr; }

    void CopyFrom(Grid* grid);
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

    bool IsAField(const Field* field) const;
    void AddField(Field* field);
    void DeleteField(Field* field);
    void ResetFields(const std::unordered_map<std::string, Field*>* fields);
    
    /**@}*/

    /**
     * @name Ghost management
     * 
     * @{
     */
    inline lid_t NGhostFront() const {
        m_assert(interp_ != nullptr, "interp cannot be null");
        return interp_->nghost_front();
    }
    inline lid_t NGhostBack() const {
        m_assert(interp_ != nullptr, "interp cannot be null");
        return interp_->nghost_back();
    }
    void GhostPull(Field* field);
    void GhostPull_Post(Field* field, const sid_t ida);
    void GhostPull_Wait(Field* field, const sid_t ida);
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

    void SetTol(const real_t refine_tol, const real_t coarsen_tol);
    void SetRecursiveAdapt(const bool recursive_adapt) { recursive_adapt_ = recursive_adapt; }

    void Refine(Field* field);
    void Coarsen(Field* field);
    void Adapt(Field* field);
    void Adapt(Field* field, SetValue* expression);
    void Adapt(std::list<Patch>* patches);

    void Adapt(void* criterion_ptr, void* interp_ptr, cback_coarsen_citerion_t coarsen_crit, cback_refine_criterion_t refine_crit, cback_interpolate_t interp);
    /**@}*/
};

#endif  // SRC_GRID_HPP_
