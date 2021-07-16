#ifndef SRC_FIELD_HPP_
#define SRC_FIELD_HPP_

#include <iostream>
#include <string>

#include "core/macros.hpp"
#include "core/types.hpp"

#define M_NFACES 6

/**
 * @brief A field contains all the information required ot access and handle a physical field
 * 
 */
class Field {
   private:
    const lda_t       lda_  = 0;          //!< indicate how many dimension is in the field [0,lda_[
    const std::string name_ = "default";  //!< the name of the field, used throughout the field management, must be unique

    bool   is_temp_      = false;   //!< indicate if a field is a temporary field and therefore will not be interpolated during the adaptation
    bool   ghost_status_ = false;   //!< indicate if the field has up-to-date ghosts or not
    bidx_t ghost_len_[2] = {0, 0};  //!< the number of ghostpoints actually up to date

    bctype_t* bctype_[M_NFACES];  //!< the boundary conditions for every direction; [X- X+ Y- Y+ Z- Z+]

   public:
    Field(std::string name, sid_t lda);
    ~Field();

    /**
     * @name get Field information
     * 
     * @{
     */
    [[nodiscard]] inline lda_t       lda() const { return lda_; }
    [[nodiscard]] inline std::string name() const { return name_; }
    [[nodiscard]] inline bool        is_temp() const { return is_temp_; }
    [[nodiscard]] inline bidx_t      get_ghost_len(const lda_t id) const { return (ghost_len_[id] * ghost_status_); }
    [[nodiscard]] inline bctype_t    bctype(const lda_t ida, const iface_t iface) const { return bctype_[iface][ida]; }
    [[nodiscard]] inline bctype_t*   bctype(const iface_t iface) const { return bctype_[iface]; }

    [[nodiscard]] inline bool ghost_status(const bidx_t ghost_len[2]) const {
        const bool is_trivial = (ghost_len[0] == 0) && (ghost_len[1] == 0);
        const bool is_large_enough = (ghost_len[0] <= ghost_len_[0]) && (ghost_len[1] <= ghost_len_[1]);
        return is_trivial || (ghost_status_ && is_large_enough) ;
    }
    /** @} */

    /**
     * @name set Field information
     * 
     * @{
     */

    // void ghost_status(const bool status) { ghost_status_ = status; }
    void ghost_len(const bidx_t* ghost_len) {
        ghost_len_[0] = ghost_len[0];
        ghost_len_[1] = ghost_len[1];
        ghost_status_ = (ghost_len[0] + ghost_len[1]) > 0;
    }
    void is_temp(const bool status) { is_temp_ = status; }
    void bctype(bctype_t* type, const iface_t iface);
    void bctype(bctype_t type, const sid_t ida, const iface_t loc);
    void bctype(bctype_t type, const sid_t ida);
    void bctype(bctype_t type);
    /** @} */
};

#endif  // SRC_FIELD_HPP_
