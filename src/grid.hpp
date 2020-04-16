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
#include "prof.hpp"

using std::list;
using std::numeric_limits;

class Grid : public ForestGrid {
   protected:
    // Field* working_callback_field_ = nullptr; /**< defines a current working field, used in the callback functions */

    map<string, Field*> fields_; /**< map of every field registered to this grid. The key is the field name (for easier access) */

    Prof*         prof_   = nullptr; /**< the profiler */
    Ghost*        ghost_  = nullptr;
    Interpolator* interp_ = nullptr;

   public:
    Grid(const lid_t ilvl, const bool isper[3], const lid_t l[3], MPI_Comm comm, Prof* prof);
    ~Grid();

    size_t LocalMemSize() const;
    size_t LocalNumDof() const;
    size_t GlobalNumDof() const;

    Interpolator* interp() { return interp_; }

    /**
     * @name Fields management
     * 
     * @{
     */
    // Field* working_callback_field() const { return working_callback_field_; }

    map<string, Field*>::const_iterator FieldBegin() const { return fields_.begin(); }
    map<string, Field*>::const_iterator FieldEnd() const { return fields_.end(); }

    lid_t NField() const { return (lid_t)(fields_.size()); }

    bool IsAField(const Field* field) const;
    void AddField(Field* field);
    void DeleteField(Field* field);

    /**@}*/

    /**
     * @name Ghost management
     * 
     * @{
     */
    void GhostPull(Field* field);
    /**@}*/

     /**
     * @name Grid adaptation
     * 
     * @{
     */
    void Refine(const sid_t delta_level);
    void Coarsen(const sid_t delta_level);
    /**@}*/

   private:
    void LoopOnGridBlock_(const bop_t op, Field* field) const;
};

#endif  // SRC_GRID_HPP_