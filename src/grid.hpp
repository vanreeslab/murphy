#ifndef SRC_GRID_HPP_
#define SRC_GRID_HPP_

#include <p8est_mesh.h>

#include <iostream>
#include <limits>
#include <map>
#include <string>

#include "field.hpp"
#include "ghost.hpp"
#include "gridblock.hpp"
#include "interpolator.hpp"
#include "murphy.hpp"
#include "patch.hpp"
#include "prof.hpp"

using std::list;
using std::numeric_limits;

/**
 * @brief implements the grid management and the related responsabilities on top of ForestGrid
 * 
 */
class Grid : public ForestGrid {
   protected:
    map<string, Field*> fields_;  //!< map of every field registered to this grid (the key is the field name, `field->name()`)

    Prof*         prof_   = nullptr;  //!< the profiler to use, may stay null
    Ghost*        ghost_  = nullptr;  //!< the ghost structure that handles one dimension of a field
    Interpolator* interp_ = nullptr;  //!< the interpolator to use for all the multilevel operations
    Interpolator* detail_ = nullptr;  //!< the interpolator to use for the detail computation

    real_t rtol_ = 1.0e-2;  //!< refinement tolerance, see @ref SetTol()
    real_t ctol_ = 1.0e-4;  //!< coarsening tolerance, see @ref SetTol()

    void* tmp_ptr_ = nullptr;  //!< temporary pointer, needed by the adaptation of the grid

   public:
   
    explicit Grid();
    Grid(const lid_t ilvl, const bool isper[3], const lid_t l[3], MPI_Comm comm, Prof* prof);
    ~Grid();

    size_t LocalMemSize() const;
    size_t LocalNumDof() const;
    size_t GlobalNumDof() const;

    Interpolator* interp() { return interp_; }
    Interpolator* detail() { return detail_; }

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

    map<string, Field*>::const_iterator FieldBegin() const { return fields_.begin(); }
    map<string, Field*>::const_iterator FieldEnd() const { return fields_.end(); }

    bool IsAField(const Field* field) const;
    void AddField(Field* field);
    void DeleteField(Field* field);
    void ResetFields(const map<string, Field*>* fields);
    /**@}*/

    /**
     * @name Ghost management
     * 
     * @{
     */
    void GhostPull(Field* field);
    void GhostPullSend(Field* field, const sid_t ida);
    void GhostPullRecv(Field* field, const sid_t ida);
    void GhostPullFill(Field* field, const sid_t ida);
    /**@}*/

    /**
     * @name Grid adaptation
     * 
     * @{
     */
    real_t rtol() const { return rtol_; }
    real_t ctol() const { return ctol_; }
    void* tmp_ptr() const { return tmp_ptr_; }

    void SetTol(const real_t refine_tol, const real_t coarsen_tol);
    void Refine(const sid_t delta_level);
    void Coarsen(const sid_t delta_level);
    void Adapt(Field* field);
    void Adapt(list<Patch>* patches);
    /**@}*/

   private:
    void LoopOnGridBlock_(const gbop_t op, Field* field) const;
};

#endif  // SRC_GRID_HPP_
