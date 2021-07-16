// #ifndef SRC_MGFAMILY_HPP_
// #define SRC_MGFAMILY_HPP_

// #include "gridblock.hpp"
// #include "wavelet.hpp"
// #include "core/macros.hpp"
// #include "core/types.hpp"

// /**
//  * @brief A familly is an association of one parent and 8 children
//  * 
//  * The familly takes full responsability of the created parents, as nobody else knows
//  * that they have been created by the coarsening of one grid
//  * 
//  */
// class MGFamily {
//    protected:
//     lid_t       parent_count_ = 0;
//     GridBlock** parents_      = nullptr;  //!< refer to my parent, maybe nullptr
//     GridBlock** children_     = nullptr;  //!< refer to my children, maybe nullptr
//    public:
//     MGFamily(const lid_t num_children);
//     ~MGFamily();

//     void AddMembers(GridBlock* parent, GridBlock* children[P8EST_CHILDREN]);

//     inline lid_t parent_count() const { return parent_count_; }

//     void ToParents(Field* field_src, Field* field_trg, Field* field_cst, const real_t alpha, Wavelet* interp);
//     void ToChildren(Field* field_src, Field* field_trg, Field* field_cst, const real_t alpha, Wavelet* interp);
// };

// /**
//  * @brief pointer to an member function of the class @ref MGFamily
//  */
// // using mgfop_t = void (MGFamily::*)();

// #endif  // SRC_MGFAMILY_HPP_
