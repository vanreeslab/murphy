#include "field.hpp"

/**
 * @brief Construct a new Field given a unique name and a dimension
 * 
 * @param name the name has to be unique
 * @param lda the dimension of the field
 */
Field::Field(string name, sid_t lda) {
    name_         = name;
    lda_          = lda;
    ghost_status_ = false;

    // create an empty BC array
    for (sid_t id = 0; id < 6; id++) {
        bctype_[id] = reinterpret_cast<bctype_t*>(m_calloc(sizeof(bctype_t) * lda));
        for(sid_t ida=0; ida<lda; ida++){
            bctype_[id][ida] = M_BC_NONE;
        }
    }
}

/**
 * @brief Destroy the Field
 * 
 */
Field::~Field() {
    // create an empty BC array
    for (sid_t id = 0; id < 6; id++) {
        m_free(bctype_[id]);
    }
}

/**
 * @brief set the same boundary condition in every direction, for a given dimension of the field
 * 
 * @warning the boundary condition will not be taken into account if the direction is periodic
 * 
 * @param type the boundary condition to set in every direction
 * @param ida the dimension of the field
 */
void Field::bctype(bctype_t type, const sid_t ida) {
    for (int iloc = 0; iloc < 6; iloc++) {
        bctype_[iloc][ida] = type;
    }
}

/**
 * @brief set the same boundary condition in every direction for every dimension of the field
 * 
 * @warning the boundary condition will not be taken into account if the direction is periodic
 * 
 * @param type the boundary condition to set in every direction for every dimension
 */
void Field::bctype(bctype_t type) {
    for (int ida = 0; ida < lda_; ida++) {
        for (int iloc = 0; iloc < 6; iloc++) {
            bctype_[iloc][ida] = type;
        }
    }
}

/**
 * @brief set a boundary condition in a given direction for a given dimension of a field
 * 
 * @warning the boundary condition will not be taken into account if the direction is periodic
 * 
 * @param type the boundary condition to set
 * @param ida the dimension of the field that will get the boundary condition
 * @param loc the placement of the boundary condition (X- = 0, X+ = 1, Y- = 2, Y+ = 3, Z- = 4, Z+ = 5 )
 */
void Field::bctype(bctype_t type, const sid_t ida, const iface_t loc) {
    bctype_[loc][ida] = type;
}

/**
 * @brief replaces the boundary condition pointer but the one given
 * 
 * @param type the new boundary condition pointer
 * @param iface the face you wish to replace
 */
void Field::bctype(bctype_t* type, const iface_t iface) {
    bctype_[iface] = type;
}
