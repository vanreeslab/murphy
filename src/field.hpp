#ifndef SRC_FIELD_HPP_
#define SRC_FIELD_HPP_

#include <iostream>
#include <string>

#include "murphy.hpp"

using std::string;
using std::to_string;

class Field {
   public:
    sid_t  lda_          ;
    string name_         ;
    bool   ghost_status_ ;

    Field(string name, sid_t lda) {
        name_ = name;
        lda_  = lda;
        ghost_status_ = false;
    }

    string name() const;
    sid_t  lda() const;

    void SetGhostStatus(bool status) { ghost_status_ = status; }
};

class SubField : public Field {
};

#endif  // SRC_FIELD_HPP_
