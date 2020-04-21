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

template <int npoint>
class Boundary {
   protected:
    const real_t face_sign_[6][3]  = {{-1.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, -1.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, -1.0}, {0.0, 0.0, 1.0}};
    const lid_t  face_start_[6][3] = {{0, 0, 0}, {M_N - 1, 0, 0}, {0, 0, 0}, {0, M_N - 1, 0}, {0, 0, 0}, {0, 0, M_N - 1}};

   public:
    virtual real_t Stencil_(const real_t* f, const real_t x, const real_t h, const real_t offset, const real_t normal, const real_t flux) = 0;

    virtual void operator()(const real_t boundary_condition, const real_t hgrid[3], PhysBlock* block, real_p data) {
        const lid_t*  fstart = face_start_[block->iface()];
        const real_t* fsign  = face_sign_[block->iface()];
        // get the face direction and a boolean array with it
        const sid_t dir       = block->iface() / 2;
        const bool  isphys[3] = {dir == 0, dir == 1, dir == 2};
        // get the offset in memory, i.e. the position of the first point (0,0,0)
        real_t offset[3];
        m_pos_relative(offset, 0, 0, 0, hgrid);

        // shift the data to the correct spot
        real_p ldata = data + m_idx(fstart[0], fstart[1], fstart[2]);

        // let's goooo
        for (lid_t i2 = block->start(2); i2 < block->end(2); i2++) {
            for (lid_t i1 = block->start(1); i1 < block->end(1); i1++) {
                for (lid_t i0 = block->start(0); i0 < block->end(0); i0++) {
                    // we need three interpolations points in the face direction
                    real_t f[npoint];
                    for (sid_t ip = 0; ip < npoint; ip++) {
                        // if the direction is not physics, just consider the current ID,
                        // if the direction is physics, takes the first, second and third point INSIDE the block (0 = first inside)
                        f[ip] = ldata[m_idx((!isphys[0]) * i0 - ip * fsign[0] * isphys[0],
                                            (!isphys[1]) * i1 - ip * fsign[1] * isphys[1],
                                            (!isphys[2]) * i2 - ip * fsign[2] * isphys[2])];
                    }
                    // get the ghost point position
                    const real_t x = (i0 * isphys[0] + i1 * isphys[1] + i2 * isphys[2]) * hgrid[dir];
                    // get the ghost value
                    ldata[m_idx(i0, i1, i2)] = Stencil_(f, x, hgrid[dir], offset[dir], fsign[dir], boundary_condition);
                }
            }
        }
    }
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
class EvenBoundary_4 : public Boundary<3> {
    real_t stencil_[3];
   protected:
    real_t Stencil_(const real_t* f, const real_t x, const real_t h, const real_t offset, const real_t normal, const real_t flux) override;
};

/**
 * @brief applies an ODD bondary condition with respect to a given value at the interface
 *
 * given a value, i.e. a constant fucntion, we impose a zero additional value and a anti-symmetry.
 * This generalizes the traditional odd symmetry condition (value = 0).
 *
 * Due to the odd symmetry, every even exponent is null in the polynomial:
 * f(x) = a x^5 + b x^3 + c x = data - value
 * whith
 *      f(h/2)  - value = a  h^5/32 + b h^3/8 + c h/2
 *      f(3h/2) - value = a  243*h^4/32 + b 27*h^3/8 + c 3*h/2
 *      f(5h/2) - value = a  3125*h^4/32 + b 125*h^3/8 + c 5*h/2
 */
class OddBoundary_4 : public Boundary<3> {
    real_t stencil[3];

   protected:
    real_t Stencil_(const real_t* f, const real_t x, const real_t h, const real_t offset, const real_t normal, const real_t value) override;
};

/**
 * @brief applies an UNBOUNDED bondary condition by extrapolation
 *
 * f(x) = a x^3 + b x^2 + c x + d = data;
 */
class ExtrapBoundary_4 : public Boundary<4> {
    real_t stencil[3];

   protected:
    real_t Stencil_(const real_t* f, const real_t x, const real_t h, const real_t offset, const real_t normal, const real_t value) override;
};


/**
 * @brief applies an UNBOUNDED bondary condition by extrapolation
 *
 * f(x) = a x^3 + b x^2 + c x + d = data;
 */
class ZeroBoundary : public Boundary<0> {
    real_t stencil[3];

   protected:
    real_t Stencil_(const real_t* f, const real_t x, const real_t h, const real_t offset, const real_t normal, const real_t value) override;
};


#endif  // SRC_BOUNDARY_HPP_
