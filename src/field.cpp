#include "field.hpp"

Field::Field(string name, sid_t lda) {
    name_         = name;
    lda_          = lda;
    ghost_status_ = false;
    // ida_          = -1;
}

// /**
//  * @brief select a working component for the field
//  * 
//  * @param ida 
//  * @return sid_t 
//  */
// void Field::SelectDimension(const sid_t ida) {
//     m_assert((0 <= ida), "the selected dimension must be positive");
//     m_assert((ida < lda_), "the selected dimension must be < lda");
//     ida_ = ida;
// }

// /**
//  * @brief Release a working direction for the field
//  * 
//  * @param ida 
//  */
// void Field::ReleaseDimension(const sid_t ida) {
//     m_assert(ida_ != -1, "no dimension has been selected beforehand");
//     ida_ = -1;
// }

// /**
//  * @brief returns the first dimension to work on
//  * 
//  * if a working dimension has been selected, using SelectDimension(), return it
//  * if no working dimension is available, start to work on the first one (=0)
//  * 
//  * @return sid_t 
//  */
// sid_t Field::StartDimension() const{
//     return (ida_ == -1) ? 0 : ida_;
// }

// /**
//  * @brief returns the last dimension (+1) to work on
//  * 
//  * if a working dimension has been selected, using SelectDimension(), return it+1
//  * if no working dimension is available, return the total number of dimension (lda_)
//  * 
//  * @return sid_t 
//  */
// sid_t Field::EndDimension() const {
//     return (ida_ == -1) ? lda_ : ida_ + 1;
// }
