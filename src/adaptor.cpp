#include "adaptor.hpp"
#include "gridcallback.hpp"
#include "p8est_extended.h"

Adaptor::Adaptor(ForestGrid* grid) {
    m_begin;
    //-------------------------------------------------------------------------
    grid_ = grid;

    //-------------------------------------------------------------------------
    m_end;
}

void Adaptor::Refine(Field* field, Interpolator* interp){
    m_begin;
    //-------------------------------------------------------------------------
    interp_ = interp;

    // take the forest and refer the current object
    p8est_t* forest = grid_->forest();

    p8est_refine_ext(forest, 0, P8EST_MAXLEVEL, cback_Yes, NULL, cback_Interpolate);

    //-------------------------------------------------------------------------
    m_end;

}