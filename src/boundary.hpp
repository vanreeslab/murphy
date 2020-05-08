#ifndef SRC_BOUNDARY_HPP_
#define SRC_BOUNDARY_HPP_

#include "murphy.hpp"
#include "physblock.hpp"

static real_t face_sign[6][3]  = {{-1.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, -1.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, -1.0}, {0.0, 0.0, 1.0}};
static lid_t  face_start[6][3] = {{0, 0, 0}, {M_N - 1, 0, 0}, {0, 0, 0}, {0, M_N - 1, 0}, {0, 0, 0}, {0, 0, M_N - 1}};

/**
 * @brief implements a boundary conditions that uses points inside the block to extrapolate a ghost value
 * 
 * @tparam npoint the number of points needed within the block to be used in @ref Stencil_()
 */
template <int npoint>
class Boundary {
   public:
   /**
    * @brief has to be implemented by boundary condition classes
    * 
    * given the data values, builds a polynomial representation of the ghost points:
    * 
    * ```
    *       GHOST   |                       BLOCK                      |    GHOST
    * ------------------------------------------------------------------------------------
    *   normal = -1 |                                                  | normal = +1
    *               |                                                  |
    *        x (<0) | offset                                           |   x (>0) 
    *     <-------->|<->                                               |<-------->
    * ---o------o---|---o------o------o---   ...  ---o------o------o---|---o------o------
    *                  f[0]   f[1]   ...            ...   f[1]    f[0]
    * ```
    * 
    * @param f the field values used to extrapolate (array of size @ref npoint, where f[0] is the closest point to the interface)
    * @param x the position of the ghost point to compute, with respect ot the interface
    * @param h the grid space separating the data values
    * @param offset the offset used in the grid spacing, i.e. the position of the point (0,0,0) wrt to the block's corner (left bottom)
    * @param normal the normal of the interface: +1 means the ghost point is on the + side, -1 means the ghost point is on the - side
    * @param value the boundary condition value, evaluated at the interface
    * @return real_t the value to set to the ghost point
    */
    virtual real_t Stencil_(const real_t* f, const real_t x, const real_t h, const real_t offset, const real_t normal, const real_t value) = 0;

    /**
     * @brief Given a subblock and a boundary condition value, apply a physical BC on a face
     * 
     * note: unless other functions using SubBlocks, we need to decouple the starting point (fstart)
     * and the symmetry done on the left and on the right of it for the interpolation, hence we need both arguments
     * 
     * @param iface the index of the face on which we apply it (x-=0, x+=1, y-=2, y+=3, z-=4, z-=5)
     * @param fstart the starting index of the face, i.e. the index of the last point wich is adjacent to the face (see @ref face_start)
     * @param hgrid the local grid spacing
     * @param boundary_condition the boundary condition value
     * @param block the sub-block on which we apply it
     * @param data the data pointer that refers to the (0,0,0) location associated to the subblock.
     */
    virtual void operator()(const sid_t iface, const lid_t fstart[3], const real_t hgrid[3], const real_t boundary_condition, SubBlock* block, real_p data) {
        m_assert(iface >= 0 && iface < 6, "iface = %d is wrong", iface);
        //-------------------------------------------------------------------------
        // get the face direction and a boolean array with it
        const sid_t   dir       = iface / 2;
        const bool    isphys[3] = {dir == 0, dir == 1, dir == 2};
        const real_t* fsign     = face_sign[iface];
        // get the offset in memory, i.e. the position of the first point (0,0,0)
        real_t offset[3];
        m_pos_relative(offset, 0, 0, 0, hgrid);

        // shift the data to the correct spot
        real_p ldata = data + m_midx(fstart[0], fstart[1], fstart[2], 0, block);

        // let's goooo
        for (lid_t i2 = block->start(2); i2 < block->end(2); i2++) {
            for (lid_t i1 = block->start(1); i1 < block->end(1); i1++) {
                for (lid_t i0 = block->start(0); i0 < block->end(0); i0++) {
                    // we need three interpolations points in the face direction
                    real_t f[npoint];
                    for (sid_t ip = 0; ip < npoint; ip++) {
                        // if the direction is not physics, just consider the current ID,
                        // if the direction is physics, takes the first, second and third point INSIDE the block (0 = first inside)
                        const lid_t idx0 = (!isphys[0]) * i0 - ip * (lid_t)fsign[0] * isphys[0];
                        const lid_t idx1 = (!isphys[1]) * i1 - ip * (lid_t)fsign[1] * isphys[1];
                        const lid_t idx2 = (!isphys[2]) * i2 - ip * (lid_t)fsign[2] * isphys[2];

                        m_assert((fstart[0] + idx0) >= (-block->gs()) && (fstart[0] + idx0) < (block->stride() - block->gs()), "index 0 is wrong: %d = %d - %d * %d * %d with gs=%d and strid=%d", fstart[0] + idx0, (!isphys[0]) * i0, ip, (int)fsign[0], (int)isphys[0], block->gs(), block->stride());
                        m_assert((fstart[1] + idx1) >= (-block->gs()) && (fstart[1] + idx1) < (block->stride() - block->gs()), "index 1 is wrong: %d = %d - %d * %d * %d with gs=%d and strid=%d", fstart[1] + idx1, (!isphys[1]) * i1, ip, (int)fsign[1], (int)isphys[1], block->gs(), block->stride());
                        m_assert((fstart[2] + idx2) >= (-block->gs()) && (fstart[2] + idx2) < (block->stride() - block->gs()), "index 2 is wrong: %d = %d - %d * %d * %d with gs=%d and strid=%d", fstart[2] + idx2, (!isphys[2]) * i2, ip, (int)fsign[2], (int)isphys[2], block->gs(), block->stride());

                        f[ip] = ldata[m_midx(idx0, idx1, idx2, 0, block)];
                    }
                    // get the ghost point position
                    real_t pos[3];
                    // if we have a negative normal, simply the position of the indexes
                    // if we have a positive normal, we need to substract 1
                    m_pos_relative(pos, i0 - (fsign[0] == +1), i1 - (fsign[1] == +1), i2 - (fsign[2] == +1), hgrid);
                    // get the ghost value
                    ldata[m_midx(i0, i1, i2, 0, block)] = Stencil_(f, pos[dir], hgrid[dir], offset[dir], fsign[dir], boundary_condition);
                }
            }
        }
        fflush(stdout);
    }
    //-------------------------------------------------------------------------
};

/**
 * @brief applies an EVEN bondary condition with respect to a given flux at the interface
 * 
 * given a flux, i.e. a slope, we impose a zero additional flux and a symmetry.
 * This generalizes the traditional even symmetry condition (flux = 0).
 * 
 * Due to the even symmetry, every odd exponent is null in the polynomial:
 * ```
 * f(x) = a x^4 + b x^2 + c = data - flux * x
 * ```
 * whith   
 * ```
 *      f(h/2)  - flux * h/2  = a  h^4/16 + b h^2/4 + c
 *      f(3h/2) - flux * 3h/2 = a  81*h^4/16 + b 9*h^2/4 + c
 *      f(5h/2) - flux * 5h/2 = a  625*h^4/16 + b 25*h^2/4 + c
 * ```
 */
class EvenBoundary_4 : public Boundary<3> {
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
 * ```
 * f(x) = a x^5 + b x^3 + c x = data - value
 * ```
 * whith
 * ```
 *      f(h/2)  - value = a  h^5/32 + b h^3/8 + c h/2
 *      f(3h/2) - value = a  243*h^4/32 + b 27*h^3/8 + c 3*h/2
 *      f(5h/2) - value = a  3125*h^4/32 + b 125*h^3/8 + c 5*h/2
 * ```
 */
class OddBoundary_4 : public Boundary<3> {
   protected:
    real_t Stencil_(const real_t* f, const real_t x, const real_t h, const real_t offset, const real_t normal, const real_t value) override;
};

/**
 * @brief applies an extrapolation boundary condition O(h^4) using a polynomial a x^3 + b x^2 + c x + d
 *
 */
class ExtrapBoundary_4 : public Boundary<4> {
   protected:
    real_t Stencil_(const real_t* f, const real_t x, const real_t h, const real_t offset, const real_t normal, const real_t value) override;
};
/**
 * @brief applies an extrapolation boundary condition O(h^3) using a polynomial a x^2 + b x + c
 */
class ExtrapBoundary_3 : public Boundary<3> {
   protected:
    real_t Stencil_(const real_t* f, const real_t x, const real_t h, const real_t offset, const real_t normal, const real_t value) override;
};
/**
 * @brief applies an extrapolation boundary condition O(h^5) using a polynomial a x^4 + b x^3 + c x^2 + d x + e
 */
class ExtrapBoundary_5 : public Boundary<5> {
   protected:
    real_t Stencil_(const real_t* f, const real_t x, const real_t h, const real_t offset, const real_t normal, const real_t value) override;
};

/**
 * @brief applies an UNBOUNDED bondary condition by extrapolation
 * ```
 * f(x) = a x^3 + b x^2 + c x + d = data;
 * ```
 */
class ZeroBoundary : public Boundary<0> {
   protected:
    real_t Stencil_(const real_t* f, const real_t x, const real_t h, const real_t offset, const real_t normal, const real_t value) override;
};

#endif  // SRC_BOUNDARY_HPP_