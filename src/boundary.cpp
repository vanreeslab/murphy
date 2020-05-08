#include "boundary.hpp"

/**
 * @brief implements Boundary::Stencil_() for @ref EvenBoundary_4
 */
real_t EvenBoundary_4::Stencil_(const real_t* f, const real_t x, const real_t h, const real_t offset, const real_t normal, const real_t flux) {
    m_assert(offset == (0.5 * h), "this function is only for cell-centered ghosts");
    // remove the flux slope form the value
    // if x is positive, we need go in the negative and vice-versa
    real_t cor_f[3] = {f[0] - normal * flux * (0.0 + offset) * h,
                       f[1] - normal * flux * (1.0 + offset) * h,
                       f[2] - normal * flux * (2.0 + offset) * h};
    // compute the stencil coefficients
    const real_t hnorm = -h * normal;
    const real_t a     = (2.0 * cor_f[0] - 3.0 * cor_f[1] + 1.0 * cor_f[2]) / (24.0 * pow(hnorm, 4));
    const real_t b     = (-34.0 * cor_f[0] + 39.0 * cor_f[1] - 5.0 * cor_f[2]) / (48.0 * pow(hnorm, 2));
    const real_t c     = (150.0 * cor_f[0] - 25.0 * cor_f[1] + 3.0 * cor_f[2]) / (128.0);
    // return the polynomial + the flux slope
    return a * pow(x, 4) + b * pow(x, 2) + c + flux * x;
}

/**
 * @brief implements Boundary::Stencil_() for @ref OddBoundary_4
 */
real_t OddBoundary_4::Stencil_(const real_t* f, const real_t x, const real_t h, const real_t offset, const real_t normal, const real_t value) {
    m_assert(offset == (0.5 * h), "this function is only for cell-centered ghosts");
    // correct the data entry
    real_t cor_f[3] = {f[0] - value,
                       f[1] - value,
                       f[2] - value};
    // compute the h coefficient, given the sign of the normal coefficients
    const real_t hnorm = -h * normal;
    const real_t a     = (10.0 * cor_f[0] - 5.0 * cor_f[1] + 1.0 * cor_f[2]) / (60.0 * pow(hnorm, 5));
    const real_t b     = (-34.0 * cor_f[0] + 13.0 * cor_f[1] - 1.0 * cor_f[2]) / (24.0 * pow(hnorm, 3));
    const real_t c     = (2250.0 * cor_f[0] - 125.0 * cor_f[1] + 9.0 * cor_f[2]) / (960.0 * pow(hnorm, 1));
    // return the polynomial + the value constant
    return a * pow(x, 5) + b * pow(x, 3) + c * pow(x, 1) + value;
}

/**
 * @brief implements Boundary::Stencil_() for @ref ExtrapBoundary_4
 */
real_t ExtrapBoundary_4::Stencil_(const real_t* f, const real_t x, const real_t h, const real_t offset, const real_t normal, const real_t value) {
    m_assert(offset == (0.5 * h), "this function is only for cell-centered ghosts");  // compute the h coefficient, given the sign of the normal coefficients
    const real_t hnorm = -h * normal;
    // compute the stencil coefficients
    const real_t a = (-1.0 * f[0] + 3.0 * f[1] - 3.0 * f[2] + 1.0 * f[3]) / (6.0 * pow(hnorm, 3));
    const real_t b = (5.0 * f[0] - 13.0 * f[1] + 11.0 * f[2] - 3.0 * f[3]) / (4.0 * pow(hnorm, 2));
    const real_t c = (-71.0 * f[0] + 141.0 * f[1] - 93.0 * f[2] + 23.0 * f[3]) / (24.0 * hnorm);
    const real_t d = (35.0 * f[0] - 35.0 * f[1] + 21.0 * f[2] - 5 * f[3]) / (16.0);
    // return the polynomial + the flux slope
    return a * pow(x, 3) + b * pow(x, 2) + c * pow(x, 1) + d;
}

/**
 * @brief implements Boundary::Stencil_() for @ref ExtrapBoundary_4
 */
real_t ExtrapBoundary_5::Stencil_(const real_t* f, const real_t x, const real_t h, const real_t offset, const real_t normal, const real_t value) {
    m_assert(offset == (0.5 * h), "this function is only for cell-centered ghosts");  // compute the h coefficient, given the sign of the normal coefficients
    const real_t hnorm = -h * normal;
    // compute the stencil coefficients
    const real_t a = (1.0 * f[0] - 4.0 * f[1] + 6.0 * f[2] - 4.0 * f[3] + 1.0 * f[4]) / (24.0 * pow(hnorm, 4));
    const real_t b = (-3.0 * f[0] + 11.0 * f[1] - 15.0 * f[2] + 9.0 * f[3] - 2.0 * f[4]) / (6.0 * pow(hnorm, 3));
    const real_t c = (103.0 * f[0] - 328.0 * f[1] + 390.0 * f[2] - 208.0 * f[3] + 43.0 * f[4]) / (48.0 * pow(hnorm, 2));
    const real_t d = (-93.0 * f[0] + 229.0 * f[1] - 225.0 * f[2] + 111.0 * f[3] - 22.0 * f[4]) / (24.0 * hnorm);
    const real_t e = (315.0 * f[0] - 420.0 * f[1] + 378.0 * f[2] - 180 * f[3] + 35.0 * f[4]) / (128.0);
    // return the polynomial + the flux slope
    return a * pow(x, 4) + b * pow(x, 3) + c * pow(x, 2)+ d * x + e;
}

/**
 * @brief implements Boundary::Stencil_() for @ref ExtrapBoundary_3
 */
real_t ExtrapBoundary_3::Stencil_(const real_t* f, const real_t x, const real_t h, const real_t offset, const real_t normal, const real_t value) {
    m_assert(offset == (0.5 * h), "this function is only for cell-centered ghosts");  // compute the h coefficient, given the sign of the normal coefficients
    const real_t hnorm = -h * normal;
    // compute the stencil coefficients
    const real_t a = (1.0 * f[0] - 2.0 * f[1] + 1.0 * f[2]) / (2.0 * pow(hnorm, 2));
    const real_t b = (-2.0 * f[0] + 3.0 * f[1] - 1.0 * f[2]) / (hnorm);
    const real_t c = (15.0 * f[0] - 10.0 * f[1] + 3.0 * f[2]) / (8.0);
    // return the polynomial + the flux slope
    return a * pow(x, 2) + b * pow(x, 1) + c;
}

/**
 * @brief implements Boundary::Stencil_() for @ref ZeroBoundary
 */
real_t ZeroBoundary::Stencil_(const real_t* f, const real_t x, const real_t h, const real_t offset, const real_t normal, const real_t value) {
    return 0.0;
}