#ifndef __INTERPOLATE_HPP
#define __INTERPOLATE_HPP

#include "murphy.hpp"
#include "p8est.h"

class Interpolate
{
public:
    // Interpolate(real_p trg,const real_p src, const dom_t sdom )=0;
protected:
    virtual void Coarsen()=0;
    virtual void Refine()=0;
    virtual void Copy()=0;

};


#endif