#ifndef SRC_ADAPTOR_HPP_
#define SRC_ADAPTOR_HPP_

#include <p8est.h>

#include "forestgrid.hpp"
#include "gridcallback.hpp"
#include "interpolator.hpp"
#include "murphy.hpp"

class Adaptor {
   protected:
    ForestGrid*   grid_;
    Interpolator* interp_;

public:
    Adaptor(ForestGrid* grid);

    void Refine(Field* field, Interpolator* interp);


};

#endif  // SRC_ADAPTOR_HPP_