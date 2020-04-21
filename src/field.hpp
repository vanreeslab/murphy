#ifndef SRC_FIELD_HPP_
#define SRC_FIELD_HPP_

#include <iostream>
#include <string>

#include "murphy.hpp"

using std::string;

class Field {
   protected:
    sid_t lda_;  //!< indicate how many dimension is in the field [0,lda_[
    // sid_t  ida_;  //!< indicate the current working dimension, if any
    string name_;
    bool   ghost_status_;

    bctype_t* bctype_[6];

   public:
    Field(string name, sid_t lda);
    ~Field();

    sid_t  lda() const { return lda_; };
    string name() const { return name_; };

    void ghost_status(bool status) { ghost_status_ = status; }
    inline bool ghost_status() { return ghost_status_ ; }

    bctype_t bctype(const sid_t ida, const sid_t iface) { return bctype_[iface][ida]; }
    void     bctype(bctype_t type, const sid_t ida, const sid_t loc) { bctype_[loc][ida] = type; }
    void     bctype(bctype_t type, const sid_t ida);
    void     bctype(bctype_t type);
};

#endif  // SRC_FIELD_HPP_
