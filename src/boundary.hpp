#ifndef SRC_BOUNDARY_HPP_
#define SRC_BOUNDARY_HPP_

#include "murphy.hpp"
#include "physblock.hpp"

// /**
//  * @brief possible positions of a Boundary condition
//  * 
//  */
// typedef enum bcloc_t {
//     M_X_MINUS, /**< X - */
//     M_X_PLUS, /**< X + */
//     M_Y_MINUS, /**< Y - */
//     M_Y_PLUS, /**< Y + */
//     M_Z_MINUS, /**< Z - */
//     M_Z_PLUS  /**< Z + */
// } bcloc_t;

class Boundary {
   public:
    virtual void   operator()(const real_t boundary_condition, const real_t hgrid[3], PhysBlock* block, real_p data);
    virtual real_t Stencil_(const real_t* f, const real_t x, const real_t h, const real_t offset, const real_t normal, const real_t flux) = 0;
};

/**
 * @brief applies an EVEN bondary condition with respect to a given flux at the interface
 * 
 * given a flux, i.e. a slope, we impose a zero additional flux and a symmetry.
 * This generalizes the traditional even symmetry condition (flux = 0).
 * 
 * Due to the even symmetry, every odd exponent is null in the polynomial:
 * f(x) = a x^4 + b x^2 + c = data - flux * x
 * whith   
 *      f(h/2)  - flux * h/2  = a  h^4/16 + b h^2/4 + c
 *      f(3h/2) - flux * 3h/2 = a  81*h^4/16 + b 9*h^2/4 + c
 *      f(5h/2) - flux * 5h/2 = a  625*h^4/16 + b 25*h^2/4 + c
 */
class EvenBoundary_4 : public Boundary {
    real_t stencil_[3];
   protected:
    real_t Stencil_(const real_t* f, const real_t x, const real_t h, const real_t offset, const real_t normal, const real_t flux) override;
};

// /**
//  * @brief applies an ODD bondary condition with respect to a given value at the interface
//  *
//  * given a value, i.e. a constant fucntion, we impose a zero additional value and a anti-symmetry.
//  * This generalizes the traditional odd symmetry condition (value = 0).
//  *
//  * Due to the even symmetry, every odd exponent is null in the polynomial:
//  * f(x) = a x^5 + b x^3 + c x = data - value
//  * whith
//  *      f(h/2)  - flux * h/2  = a  h^4/16 + b h^2/4 + c
//  *      f(3h/2) - flux * 3h/2 = a  81*h^4/16 + b 9*h^2/4 + c
//  *      f(5h/2) - flux * 5h/2 = a  625*h^4/16 + b 25*h^2/4 + c
//  */
// class OddBoundary_4 : public Boundary {
//     real_t stencil[3];

//    public:
//     virtual void operator()(SubBlock* block,const real_t bc_value) override;

//    protected:
//     real_t StencilZeroFlux_(const real_t* d, const real_t x, const real_t h);
// };

#endif  // SRC_BOUNDARY_HPP_
