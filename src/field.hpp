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

   public:
    Field(string name, sid_t lda);

    sid_t  lda() const { return lda_; };
    string name() const { return name_; };

    void ghost_status(bool status) { ghost_status_ = status; }
    inline bool ghost_status() { return ghost_status_ ; }

    // void SelectDimension(const sid_t ida);
    // void ReleaseDimension(const sid_t ida);

    // sid_t StartDimension() const;
    // sid_t EndDimension() const;
};

#endif  // SRC_FIELD_HPP_
