#ifndef SRC_FIELD_HPP_
#define SRC_FIELD_HPP_

#include <iostream>
#include <string>

#include "murphy.hpp"

using std::string;

/**
 * @brief Contains all the required information about a physical field
 * 
 */
class Field {
   protected:
    sid_t     lda_;           //!< indicate how many dimension is in the field [0,lda_[
    bool      ghost_status_;  //!< indicate if the field has up-to-date ghosts or not
    string    name_;          //!< the name of the field, used throughout the field management, must be unique
    bctype_t* bctype_[6];     //!< the boundary conditions for every direction; [X- X+ Y- Y+ Z- Z+]

   public:
    Field(string name, sid_t lda);
    ~Field();

    /**
     * @name get Field information
     * 
     * @{
     */
    inline sid_t     lda() const { return lda_; }
    inline string    name() const { return name_; }
    inline bool      ghost_status() const { return ghost_status_; }
    inline bctype_t  bctype(const sid_t ida, const sid_t iface) const { return bctype_[iface][ida]; }
    inline bctype_t* bctype(const sid_t iface) const { return bctype_[iface]; }
    /** @} */

    /**
     * @name set Field information
     * 
     * @{
     */

    /**
     * @brief sets the ghost information
     * 
     * @param status the ghost information status to set
     */
    void ghost_status(bool status) { ghost_status_ = status; }

    void bctype(bctype_t* type, const iface_t iface);
    void bctype(bctype_t type, const sid_t ida, const iface_t loc);
    void bctype(bctype_t type, const sid_t ida);
    void bctype(bctype_t type);

    /** @} */
};

#endif  // SRC_FIELD_HPP_
