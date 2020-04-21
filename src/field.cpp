#include "field.hpp"

Field::Field(string name, sid_t lda) {
    name_         = name;
    lda_          = lda;
    ghost_status_ = false;

    // create an empty BC array
    for (sid_t id = 0; id < 6; id++) {
        bctype_[id] = (bctype_t*)m_calloc(sizeof(bctype_t) * lda);
    }
}

Field::~Field() {
    // create an empty BC array
    for (sid_t id = 0; id < 6; id++) {
        m_free(bctype_[id]);
    }
}

void Field::bctype(bctype_t type, const sid_t ida) {
    for (int iloc = 0; iloc < 6; iloc++) {
        bctype_[iloc][ida] = type;
    }
}
void Field::bctype(bctype_t type) {
    for (int ida = 0; ida < lda_; ida++) {
        for (int iloc = 0; iloc < 6; iloc++) {
            bctype_[iloc][ida] = type;
        }
    }
}