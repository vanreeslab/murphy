// #include "boundary.hpp"

// /**
//  * @brief applies an EVEN bondary condition with respect to a given flux at the interface
//  * 
//  * given a flux, i.e. a slope, we impose a zero additional flux and a symmetry.
//  * This generalizes the traditional even symmetry condition (flux = 0).
//  * 
//  * Due to the even symmetry, every odd exponent is null in the polynomial:
//  * ```
//  * f(x) = a x^4 + b x^2 + c = data - flux * x
//  * ```
//  * if    
//  * ```
//  *      f(h)  - flux * h  = a  h^4/16 + b h^2/4 + c
//  *      f(2h) - flux * 2h = a  81*h^4/16 + b 9*h^2/4 + c
//  *      f(3h) - flux * 3h = a  625*h^4/16 + b 25*h^2/4 + c
//  * ```
//  */
// // real_t EvenBoundary_4::Stencil_(const real_t* f, const real_t x, const real_t h, const real_t normal, const real_t flux) {
// //     //-------------------------------------------------------------------------
// //     // if x is positive, we need go in the negative and vice-versa
// //     real_t cor_f[3] = {f[0] - normal * flux * (1.0) * h,
// //                        f[1] - normal * flux * (2.0) * h,
// //                        f[2] - normal * flux * (3.0) * h};
// //     // compute the stencil coefficients
// //     const real_t hnorm = -h * normal;
// //     // const real_t a     = (2.0 * cor_f[0] - 3.0 * cor_f[1] + 1.0 * cor_f[2]) / (24.0 * pow(hnorm, 4));
// //     // const real_t b     = (-34.0 * cor_f[0] + 39.0 * cor_f[1] - 5.0 * cor_f[2]) / (48.0 * pow(hnorm, 2));
// //     // const real_t c     = (150.0 * cor_f[0] - 25.0 * cor_f[1] + 3.0 * cor_f[2]) / (128.0);
// //     const real_t a = (2.0 * cor_f[0] - 3.0 * cor_f[1] + 1.0 * cor_f[2]) / (384.0 * pow(hnorm, 4));
// //     const real_t b = (-34.0 * cor_f[0] + 39.0 * cor_f[1] - 5.0 * cor_f[2]) / (192.0 * pow(hnorm, 2));
// //     const real_t c = (150.0 * cor_f[0] - 25.0 * cor_f[1] + 3.0 * cor_f[2]) / (128.0);
// //     // return the polynomial + the flux slope (x already has the correct sign)
// //     return a * pow(x, 4) + b * pow(x, 2) + c + flux * x;
// //     //-------------------------------------------------------------------------
// // }
// real_t EvenBoundary_4::Stencil_(real_t* f,  real_t* xf, const real_t xp, const real_t bc_value) {
//     return 0.0;
// }

// /**
//  * @brief applies an ODD bondary condition with respect to a given value at the interface
//  *
//  * given a value we impose a zero additional value and a anti-symmetry.
//  * This generalizes the traditional odd symmetry condition (value = 0).
//  *
//  * Due to the odd symmetry, every even exponent is null in the polynomial:
//  * ```
//  * f(x) = a x^5 + b x^3 + c x = data - value
//  * ```
//  * whith
//  * ```
//  *      f(h)  - value = a     h^5  + b    h^3 + c * h
//  *      f(2h) - value = a  32*h^5  + b  8*h^3 + c * 3h
//  *      f(3h) - value = a  243*h^5 + b 27*h^3 + c * 5h
//  * ```
//  */
// // real_t OddBoundary_4::Stencil_(const real_t* f, const real_t x, const real_t h, const real_t normal, const real_t value) {
// //     //-------------------------------------------------------------------------
// //     // correct the data entry
// //     real_t cor_f[3] = {f[0] - value,
// //                        f[1] - value,
// //                        f[2] - value};
// //     // compute the h coefficient, given the sign of the normal coefficients
// //     const real_t hnorm = -h * normal;
// //     const real_t a     = (10.0 * cor_f[0] - 5.0 * cor_f[1] + 1.0 * cor_f[2]) / (1920.0 * pow(hnorm, 5));
// //     const real_t b     = (-34.0 * cor_f[0] + 13.0 * cor_f[1] - 1.0 * cor_f[2]) / (192.0 * pow(hnorm, 3));
// //     const real_t c     = (2250.0 * cor_f[0] - 125.0 * cor_f[1] + 9.0 * cor_f[2]) / (1920.0 * pow(hnorm, 1));
// //     // return the polynomial + the value constant
// //     return a * pow(x, 5) + b * pow(x, 3) + c * pow(x, 1) + value;
// //     //-------------------------------------------------------------------------
// // }
// // real_t OddBoundary_4::Stencil_(real_t* f,  real_t* xf, const real_t xp, const real_t bc_value) {
// //     return 0.0;
// // }

// /**
//  * @brief applies an UNBOUNDED bondary condition by extrapolation
//  * ```
//  * f(x) = b x^2 + c x + d = data;
//  * ```
//  */
// // real_t ExtrapBoundary_3::Stencil_(real_t* f, const real_t x, const real_t h, const real_t normal, const real_t value) {
// real_t ExtrapBoundary_3::Stencil_(real_t* f,  real_t* xf, const real_t xp, const real_t bc_value){
//     //-------------------------------------------------------------------------
//     // const real_t hnorm = -h * normal;
//     // do the first round
//     // first round:
//     // p12 = (P1 * (x-x2) + (x-x1) * P2)/(x1-x2)
//     f[0] = (f[0] * (xp - xf[1]) + (xp - xf[0]) * f[1])/(xf[0]-xf[1]);
//     // p12 = (P1 * (x-x2) + (x-x1) * P2)/(x1-x2)
//     f[1] = (f[1] * (xp - xf[2]) + (xp - xf[1]) * f[2])/(xf[1]-xf[2]);
    
//     // return p012 = (P01 * (x-x2) + (x-x0) * P12)/(x0-x2)
//     return (f[0] * (xp - xf[2]) + (xp - xf[0]) * f[1])/(xf[0]-xf[2]);
//     // compute the stencil coefficients
//     // const real_t a = (1.0 * f[0] - 2.0 * f[1] + 1.0 * f[2]) / (2.0 * pow(hnorm, 2));
//     // const real_t b = (-5.0 * f[0] + 8.0 * f[1] - 3.0 * f[2]) / (2.0 * hnorm);
//     // const real_t c = (3.0 * f[0] - 3.0 * f[1] + 1.0 * f[2]);
//     // m_log("the coefs are a=%e b=%e c=%e",a,b,c);
//     // // return the polynomial + the flux slope
//     // return a * pow(x, 2) + b * pow(x, 1) + c;
//     //-------------------------------------------------------------------------
// }

// /**
//  * @brief applies an UNBOUNDED bondary condition by extrapolation
//  * ```
//  * f(x) = a x^3 + b x^2 + c x + d = data;
//  * ```
//  */
// real_t ExtrapBoundary_4::Stencil_(real_t* f,  real_t* xf, const real_t xp, const real_t bc_value) {
//     return 0.0;
// }
// // real_t ExtrapBoundary_4::Stencil_(const real_t* f, const real_t x, const real_t xp, const real_t h, const real_t normal, const real_t value) {
// //     //-------------------------------------------------------------------------
// //     const real_t hnorm = -h * normal;
// //     // compute the stencil coefficients
// //     const real_t a = (-1.0 * f[0] + 3.0 * f[1] - 3.0 * f[2] + 1.0 * f[3]) / (6.0 * pow(hnorm, 3));
// //     const real_t b = (3.0 * f[0] - 8.0 * f[1] + 7.0 * f[2] - 2.0 * f[3]) / (2.0 * pow(hnorm, 2));
// //     const real_t c = (-26.0 * f[0] + 57.0 * f[1] - 42.0 * f[2] + 11.0 * f[3]) / (6.0 * hnorm);
// //     const real_t d = (4.0 * f[0] - 6.0 * f[1] + 4.0 * f[2] - 1.0 * f[3]);
// //     // return the polynomial + the flux slope
// //     m_log("the coefs are a=%e b=%e c=%e, d=%e",a,b,c,d);
// //     return a * pow(x, 3) + b * pow(x, 2) + c * pow(x, 1) + d;
// //     //-------------------------------------------------------------------------
// // }

// /**
//  * @brief applies an UNBOUNDED bondary condition by extrapolation
//  * ```
//  * f(x) = a x^4 + b x^3 + c x^2 + d^x + e = data;
//  * ```
//  */
// real_t ExtrapBoundary_5::Stencil_(real_t* f,  real_t* xf, const real_t xp, const real_t bc_value) {
//     return 0.0;
// }
// // real_t ExtrapBoundary_5::Stencil_(const real_t* f, const real_t x, const real_t h, const real_t normal, const real_t value) {
// //     //-------------------------------------------------------------------------
// //     const real_t hnorm = -h * normal;
// //     // compute the stencil coefficients
// //     const real_t a = (1.0 * f[0] - 4.0 * f[1] + 6.0 * f[2] - 4.0 * f[3] + 1.0 * f[4]) / (24.0 * pow(hnorm, 4));
// //     const real_t b = (-7.0 * f[0] + 26.0 * f[1] - 36.0 * f[2] + 22.0 * f[3] - 5.0 * f[4]) / (12.0 * pow(hnorm, 3));
// //     const real_t c = (71.0 * f[0] - 236.0 * f[1] + 294.0 * f[2] - 164.0 * f[3] + 35.0 * f[4]) / (24.0 * pow(hnorm, 2));
// //     const real_t d = (-77.0 * f[0] + 214.0 * f[1] - 234.0 * f[2] + 122.0 * f[3] - 25.0 * f[4]) / (12.0 * hnorm);
// //     const real_t e = (5.0 * f[0] - 10.0 * f[1] + 10.0 * f[2] - 5.0 * f[3] + 1.0 * f[4]);
// //     // return the polynomial + the flux slope
// //     m_log("the coefs are a=%e b=%e c=%e, d=%e, e=%e",a,b,c,d,e);
// //     return a * pow(x, 4) + b * pow(x, 3) + c * pow(x, 2) + d * x + e;
// //     //-------------------------------------------------------------------------
// // }

// /**
//  * @brief implements Boundary::Stencil_() for @ref ZeroBoundary
//  */
// // real_t ZeroBoundary::Stencil_(const real_t* f, const real_t x, const real_t h, const real_t normal, const real_t value) {
// //     //-------------------------------------------------------------------------
// //     return 0.0;
// //     //-------------------------------------------------------------------------
// // }
// real_t ZeroBoundary::Stencil_(real_t* f,  real_t* xf, const real_t xp, const real_t bc_value) {
//     return 0.0;
// }
