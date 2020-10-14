#ifndef SRC_BOUNDARY_HPP_
#define SRC_BOUNDARY_HPP_

#include "murphy.hpp"
#include "physblock.hpp"

static real_t face_sign[6][3] = {{-1.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, -1.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, -1.0}, {0.0, 0.0, 1.0}};
/**
 * @brief the localization of the interface
 */
static lid_t face_start[6][3] = {{0, 0, 0}, {M_N, 0, 0}, {0, 0, 0}, {0, M_N, 0}, {0, 0, 0}, {0, 0, M_N}};

/**
 * @brief implements a boundary conditions that uses points inside the block to extrapolate a ghost value
 * 
 * @tparam npoint the number of points needed within the block to be used in @ref Stencil_()
 */
template <lda_t npoint>
class Boundary {
   public:
    /**
    * @brief has to be implemented by boundary condition classes
    * 
    * given the data values, builds a polynomial representation of the ghost points:
    * 
    * ```
    *          GHOST   |                       BLOCK                      |    GHOST
    * -----------------|--------------------------------------------------|---------------
    *      normal = -1 |                                                  | normal = +1
    *                  |                                                  |
    *           x (<0) |                                                  |   x (>0) 
    *     <----------->|                                                  |<-------->
    * ---o------o------o------o------o---   ...       -----o------o-------o------o------
    *                       f[0]   f[1]     ...         f[1]    f[0]
    * ```
    * 
    * @param f the field values used to extrapolate (array of size @ref npoint, where f[0] is the closest point to the interface)
    * @param x the position of the ghost point to compute, with respect ot the interface
    * @param h the grid space separating the data values
    * @param normal the normal of the interface: +1 means the ghost point is on the + side, -1 means the ghost point is on the - side
    * @param value the boundary condition value, evaluated at the interface
    * @return real_t the value to set to the ghost point
    */
    // virtual real_t Stencil_(const real_t* f, const real_t x, const real_t h, const real_t normal, const real_t value) = 0;
    virtual inline real_t Stencil_(real_t* f, real_t* xf, const real_t xp, const real_t bc_value) = 0;

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
    virtual void operator()(const sid_t iface, const lid_t fstart[3], const real_t hgrid[3], const real_t boundary_condition, SubBlock* block, data_ptr data) {
        m_assert(iface >= 0 && iface < 6, "iface = %d is wrong", iface);
        //-------------------------------------------------------------------------
        // get the face direction and a boolean array with it
        const sid_t   dir       = iface / 2;
        const bool    isphys[3] = {dir == 0, dir == 1, dir == 2};
        const real_t* fsign     = face_sign[iface];
        m_assert((fsign[0] + fsign[1] + fsign[2] == -1) || (fsign[0] + fsign[1] + fsign[2] == +1), "only 1 component of face_sign must be non null: %f %f %f ", fsign[0], fsign[1], fsign[2]);

        // shift the data to compute the region on the left or on the right of fstart
        const lid_t start[3] = {block->start(0) - fstart[0], block->start(1) - fstart[1], block->start(2) - fstart[2]};
        const lid_t end[3]   = {block->end(0) - fstart[0], block->end(1) - fstart[1], block->end(2) - fstart[2]};

        // move the data around fstart
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        data_ptr ldata = data + m_midx(fstart[0], fstart[1], fstart[2], 0, block);
        // let's goooo
        for (lid_t i2 = start[2]; i2 < end[2]; i2++) {
            for (lid_t i1 = start[1]; i1 < end[1]; i1++) {
                for (lid_t i0 = start[0]; i0 < end[0]; i0++) {
                    // we need three interpolations points in the face direction
                    real_t f[npoint];
                    real_t xf[npoint];
                    for (sid_t ip = 0; ip < npoint; ip++) {
                        // if the direction is not physics, just consider the current ID,
                        // if the direction is physics, takes the first, second and third point INSIDE the block (0 = first inside)
                        // these are local increments on the source which is already in position fstart
                        const lid_t idx0 = (!isphys[0]) ? i0 : (-(ip + 1) * fsign[0]);
                        const lid_t idx1 = (!isphys[1]) ? i1 : (-(ip + 1) * fsign[1]);
                        const lid_t idx2 = (!isphys[2]) ? i2 : (-(ip + 1) * fsign[2]);

                        // check that we stay at the correct sport for everybody
                        m_assert((fstart[0] + idx0) >= (-block->gs()) && (fstart[0] + idx0) < (block->stride()), "index 0 is wrong: %d with gs = %d and stride = %d", fstart[0] + idx0, block->gs(), block->stride());
                        m_assert((fstart[1] + idx1) >= (-block->gs()) && (fstart[1] + idx1) < (block->stride()), "index 1 is wrong: %d with gs = %d and stride = %d", fstart[1] + idx1, block->gs(), block->stride());
                        m_assert((fstart[2] + idx2) >= (-block->gs()) && (fstart[2] + idx2) < (block->stride()), "index 2 is wrong: %d with gs = %d and stride = %d", fstart[2] + idx2, block->gs(), block->stride());
                        // check the specific isphys direction
                        m_assert((((fstart[0] + idx0) * isphys[0]) >= 0) && ((fstart[0] + idx0) * isphys[0]) < (block->stride() - block->gs()), "index 0 is wrong: %d with gs = %d and stride = %d", fstart[0] + idx0, block->gs(), block->stride());
                        m_assert((((fstart[1] + idx1) * isphys[1]) >= 0) && ((fstart[1] + idx1) * isphys[1]) < (block->stride() - block->gs()), "index 1 is wrong: %d with gs = %d and stride = %d", fstart[1] + idx1, block->gs(), block->stride());
                        m_assert((((fstart[2] + idx2) * isphys[2]) >= 0) && ((fstart[2] + idx2) * isphys[2]) < (block->stride() - block->gs()), "index 2 is wrong: %d with gs = %d and stride = %d", fstart[2] + idx2, block->gs(), block->stride());

                        // store the result
                        f[ip] = ldata[m_midx(idx0, idx1, idx2, 0, block)];
                        // get the data position
                        real_t data_pos[3];
                        m_pos_relative(data_pos, idx0, idx1, idx2, hgrid);
                        xf[ip] = data_pos[dir];
                    }
                    // get the position of the ghost
                    real_t pos[3];
                    m_pos_relative(pos, i0, i1, i2, hgrid);

                    // get the ghost value
                    ldata[m_midx(i0, i1, i2, 0, block)] = Stencil_(f, xf, pos[dir], boundary_condition);
                }
            }
        }
    }
    //-------------------------------------------------------------------------
};

class ZeroBoundary : public Boundary<0> {
   protected:
    // real_t Stencil_(const real_t* f, const real_t x, const real_t h, const real_t normal, const real_t value) override;
    real_t Stencil_(real_t* f, real_t* xf, const real_t xp, const real_t bc_value) override { return 0.0; };
};


template <lda_t npoint>
class ExtrapBoundary : public Boundary<npoint> {
   protected:
    /**
     * @brief Impletement the Neville algorithm to extrapolate the value of the polynomial with no boundary information
     * 
     * The Neville implementation constructs the Lagrange polynomial in a cheap and robust way.
     * We refer to Numerical Recipes in C book for more details.
     * 
     * Compared to the version presented p109, we always have an ns = 0 as we always "extrapolate on the 0th side of the f array"
     * 
     * @param f 
     * @param xf 
     * @param xp 
     * @param bc_value 
     * @return real_t 
     */
    inline real_t Stencil_(real_t* f, real_t* xf, const real_t xp, const real_t bc_value) override {
        //-------------------------------------------------------------------------
        real_t d[npoint];
        real_t c[npoint];
        // m = 0
        for (lda_t id = 0; id < npoint; ++id) {
            d[id] = f[id];
            c[id] = f[id];
        }
        real_t value = c[0];  // ns = 0
        // m = 1 -> npoint
        for (lda_t m = 1; m < npoint; ++m) {
            for (lda_t i = 0; i < (npoint - m); ++i) {
                const real_t ho  = xf[i] - xp;
                const real_t hp  = xf[i + m] - xp;
                const real_t den = xf[i] - xf[i + m];
                const real_t w   = (c[i + 1] - d[i]) / den;
                d[i] = hp * w;
                c[i] = ho * w;
            }
            value += c[0];  // ns = 0
        }
        // m_log("return %e",f[0]);
        return value;
        //-------------------------------------------------------------------------
    }
};

template <lda_t npoint>
class OddBoundary : public Boundary<npoint> {
   protected:
    /**
     * @brief Impletement the Neville algorithm to obtain an ODD polynomial around a value, i.e. a Dirichlet boundary condition
     * 
     * We refer to Numerical Recipes in C book for more details.
     * Compared to the version presented p109, we always have an ns = 0 as we always "extrapolate on the 0th side of the f array"
     * 
     * :warning: this is quite overshoot as a method but one of the only general enough... -> the result is trivial
     * 
     * @param f 
     * @param xf 
     * @param xp 
     * @param bc_value 
     * @return real_t 
     */
    real_t Stencil_(real_t* f, real_t* xf, const real_t xp, const real_t bc_value) override {
        //-------------------------------------------------------------------------
        constexpr lda_t len  = 2 * npoint;
        constexpr lda_t sym  = npoint - 1;
        constexpr lda_t half = npoint;

        // get the arrays
        real_t d[len];
        real_t c[len];
        real_t xe[len];

        // set the rest of the vector
        for (lda_t id = 0; id < npoint; ++id) {
            const real_t fval = f[id] - bc_value;
            // store the value
            c[half + id]  = fval;
            d[half + id]  = fval;
            xe[half + id] = xf[id];
            // anti-symmetry
            d[sym - id]  = -fval;
            c[sym - id]  = -fval;
            xe[sym - id] = -xf[id];
        }

        real_t value = bc_value + c[half];  // ns = 0
        // m = 1 -> npoint
        for (lda_t m = 1; m < len; ++m) {
            for (lda_t i = 0; i < (len - m); ++i) {
                const real_t ho  = xe[i] - xp;
                const real_t hp  = xe[i + m] - xp;
                const real_t den = xe[i] - xe[i + m];
                const real_t w   = (c[i + 1] - d[i]) / den;
                d[i]             = hp * w;
                c[i]             = ho * w;
            }
            value += c[0];  // ns = 0
        }
        // m_log("return %e",f[0]);
        return value;
        //-------------------------------------------------------------------------
    }
};

template <lda_t npoint>
class EvenBoundary : public Boundary<npoint> {
   protected:
    /**
     * @brief Impletement the Neville algorithm to obtain an EVEN polynomial around a flux, i.e. a Neuman boundary condition
     * 
     * We refer to Numerical Recipes in C book for more details.
     * Compared to the version presented p109, we always have an ns = 0 as we always "extrapolate on the 0th side of the f array"
     * 
     * :warning: this is quite overshoot as a method but one of the only general enough... -> the result is trivial
     * 
     * @param f 
     * @param xf 
     * @param xp 
     * @param bc_value 
     * @return real_t 
     */
    real_t Stencil_(real_t* f, real_t* xf, const real_t xp, const real_t bc_value) override {
        //-------------------------------------------------------------------------
        constexpr lda_t len  = 2 * npoint;
        constexpr lda_t sym  = npoint - 1;
        constexpr lda_t half = npoint;

        // get the arrays
        real_t d[len];
        real_t c[len];
        real_t xe[len];

        // set the rest of the vector
        for (lda_t id = 0; id < npoint; ++id) {
            const real_t fval = f[id] - bc_value;
            // store the value
            c[half + id]  = fval;
            d[half + id]  = fval;
            xe[half + id] = xf[id];
            // anti-symmetry
            d[sym - id]  = fval;
            c[sym - id]  = fval;
            xe[sym - id] = -xf[id];
        }

        real_t value = bc_value * xp + c[half];  // ns = 0
        // m = 1 -> npoint
        for (lda_t m = 1; m < len; ++m) {
            for (lda_t i = 0; i < (len - m); ++i) {
                const real_t ho  = xe[i] - xp;
                const real_t hp  = xe[i + m] - xp;
                const real_t den = xe[i] - xe[i + m];
                const real_t w   = (c[i + 1] - d[i]) / den;
                d[i]             = hp * w;
                c[i]             = ho * w;
            }
            value += c[0];  // ns = 0
        }
        // m_log("return %e",f[0]);
        return value;
        //-------------------------------------------------------------------------
    }
};

#endif  // SRC_BOUNDARY_HPP_
