#ifndef SRC_MGFAMILY_HPP_
#define SRC_MGFAMILY_HPP_

#include "gridblock.hpp"
#include "gridblock.hpp"
#include "murphy.hpp"

#include "interpolator.hpp"

class MGFamily {
   protected:
    lid_t       parent_count_ = 0;
    GridBlock** parents_      = nullptr;  //!< refer to my parent, maybe nullptr
    GridBlock** children_     = nullptr;  //!< refer to my children, maybe nullptr
   public:
    MGFamily(const lid_t num_children);
    ~MGFamily();

    void AddMembers(GridBlock* parent, GridBlock* children[P8EST_CHILDREN]);

    inline lid_t parent_count() const { return parent_count_; }

    void ToParents(Field* field_src, Field* field_trg, Interpolator* interp);
    void ToChildren(Field* field_src, Field* field_trg, Interpolator* interp);

    // MGFamily(GridBlock* me, GridBlock* sibling);
    // MGFamily(GridBlock* parent, GridBlock* children[P8EST_CHILDREN]);
};

/**
 * @brief pointer to an member function of the class @ref MGFamily
 */
// using mgfop_t = void (MGFamily::*)();

#endif  // SRC_MGFAMILY_HPP_
