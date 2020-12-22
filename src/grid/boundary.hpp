#ifndef SRC_GRID_BOUNDARY_HPP_
#define SRC_GRID_BOUNDARY_HPP_

#include "core/macros.hpp"
#include "core/types.hpp"
#include "core/pointers.hpp"
#include "core/forloop.hpp"
#include "grid/subblock.hpp"

/**
 * @brief implements a boundary conditions that uses points inside the block to extrapolate a ghost value
 * 
 * @tparam npoint the number of points needed within the block to be used in @ref Stencil_()
 */
template <lda_t npoint>
class Boundary {
   protected:
    static constexpr real_t face_sign[6][3] = {{-1.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, -1.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, -1.0}, {0.0, 0.0, 1.0}};

    /**
     * @brief dictates if the first point in the grid, i.e. the point [0] in the physical direction is overwritten or not.
     * 
     * if the sign is different than 0, return true
     */
    virtual inline bool OverWriteFirst(const real_t sign) const { return (sign * sign) > 0.5; }

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
    * The point `x` is overwritten or not, given the boolean returned by OverWriteFirst().
    */
    virtual inline real_t Stencil_(m_ptr<const real_t> f, m_ptr<const real_t> xf, const real_t xp, const real_t bc_value) = 0;

    /**
     * @brief Given a subblock and a boundary condition value, apply a physical BC on a face
     * 
     * note: unless other functions using SubBlocks, we need to decouple the starting point (fstart)
     * and the symmetry done on the left and on the right of it for the interpolation, hence we need both arguments
     * 
     * @param iface the index of the face on which we apply it (x-=0, x+=1, y-=2, y+=3, z-=4, z-=5)
     * @param fstart the starting index of the face, i.e. the index of the last point wich is adjacent to the face (see @ref face_start in ghost.cpp)
     * @param hgrid the local grid spacing
     * @param boundary_condition the boundary condition value
     * @param block the sub-block on which we apply it
     * @param data the data pointer that refers to the (0,0,0) location associated to the subblock.
     */
    virtual void operator()(const sid_t iface, const lid_t fstart[3], const real_t hgrid[3], const real_t boundary_condition, m_ptr<const SubBlock> block, data_ptr data) {
        m_assert(iface >= 0 && iface < 6, "iface = %d is wrong", iface);
        m_assert((npoint + 1) <= M_N, "the size of the box is not big enough to take the needed samples");
        //-------------------------------------------------------------------------
        // get the face direction and a boolean array with it
        const iface_t dir       = iface / 2;
        const bool    isphys[3] = {dir == 0, dir == 1, dir == 2};
        const real_t* fsign     = face_sign[iface];
        m_assert((fsign[0] + fsign[1] + fsign[2] == -1) || (fsign[0] + fsign[1] + fsign[2] == +1), "only 1 component of face_sign must be non null: %f %f %f ", fsign[0], fsign[1], fsign[2]);

        // shift the data to compute the region on the left or on the right of fstart
        const bidx_t start[3] = {block->start(0) - fstart[0], block->start(1) - fstart[1], block->start(2) - fstart[2]};
        const bidx_t end[3]   = {block->end(0) - fstart[0], block->end(1) - fstart[1], block->end(2) - fstart[2]};

        // move the data around fstart
        rank_t rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        const bidx_t b_stride = block->stride();
        real_t*      ldata   = data.Write(fstart[0], fstart[1], fstart[2], 0, b_stride);

        auto op = [=, &ldata](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
            // we need three interpolations points in the face direction
            real_t f[npoint];
            real_t xf[npoint];
            for (bidx_t ip = 0; ip < npoint; ip++) {
                // if the direction is not physics, just consider the current ID,
                // if the direction is physics, takes the first, second and third point INSIDE the block (0 = first inside)
                // these are local increments on the source which is already in position fstart
                const bidx_t idx0 = (!isphys[0]) ? i0 : (-(ip + OverWriteFirst(fsign[0])) * fsign[0]);
                const bidx_t idx1 = (!isphys[1]) ? i1 : (-(ip + OverWriteFirst(fsign[1])) * fsign[1]);
                const bidx_t idx2 = (!isphys[2]) ? i2 : (-(ip + OverWriteFirst(fsign[2])) * fsign[2]);

                // check that we stay at the correct sport for everybody
                m_assert((fstart[0] + idx0) >= (-block->gs()) && (fstart[0] + idx0) < (block->stride()), "index 0 is wrong: %d with gs = %d and stride = %d", fstart[0] + idx0, block->gs(), block->stride());
                m_assert((fstart[1] + idx1) >= (-block->gs()) && (fstart[1] + idx1) < (block->stride()), "index 1 is wrong: %d with gs = %d and stride = %d", fstart[1] + idx1, block->gs(), block->stride());
                m_assert((fstart[2] + idx2) >= (-block->gs()) && (fstart[2] + idx2) < (block->stride()), "index 2 is wrong: %d with gs = %d and stride = %d", fstart[2] + idx2, block->gs(), block->stride());
                // check the specific isphys direction
                m_assert((((fstart[0] + idx0) * isphys[0]) >= 0) && ((fstart[0] + idx0) * isphys[0]) < (block->stride() - block->gs()), "index 0 is wrong: %d with gs = %d and stride = %d", fstart[0] + idx0, block->gs(), block->stride());
                m_assert((((fstart[1] + idx1) * isphys[1]) >= 0) && ((fstart[1] + idx1) * isphys[1]) < (block->stride() - block->gs()), "index 1 is wrong: %d with gs = %d and stride = %d", fstart[1] + idx1, block->gs(), block->stride());
                m_assert((((fstart[2] + idx2) * isphys[2]) >= 0) && ((fstart[2] + idx2) * isphys[2]) < (block->stride() - block->gs()), "index 2 is wrong: %d with gs = %d and stride = %d", fstart[2] + idx2, block->gs(), block->stride());

                // store the result
                f[ip] = ldata[m_idx(idx0, idx1, idx2, 0, b_stride)];
                // get the data position
                real_t data_pos[3];
                m_pos_relative(data_pos, idx0, idx1, idx2, hgrid);
                xf[ip] = data_pos[dir];
            }
            // get the position of the ghost
            real_t pos[3];
            m_pos_relative(pos, i0, i1, i2, hgrid);

            // get the ghost value
            ldata[m_idx(i0, i1, i2, 0, b_stride)] = Stencil_(f, xf, pos[dir], boundary_condition);
        };

        for_loop(&op,start,end);

        // // let's goooo
        // for (lid_t i2 = start[2]; i2 < end[2]; i2++) {
        //     for (lid_t i1 = start[1]; i1 < end[1]; i1++) {
        //         for (lid_t i0 = start[0]; i0 < end[0]; i0++) {
        //             // we need three interpolations points in the face direction
        //             real_t f[npoint];
        //             real_t xf[npoint];
        //             for (sid_t ip = 0; ip < npoint; ip++) {
        //                 // if the direction is not physics, just consider the current ID,
        //                 // if the direction is physics, takes the first, second and third point INSIDE the block (0 = first inside)
        //                 // these are local increments on the source which is already in position fstart
        //                 const lid_t idx0 = (!isphys[0]) ? i0 : (-(ip + OverWriteFirst(fsign[0])) * fsign[0]);
        //                 const lid_t idx1 = (!isphys[1]) ? i1 : (-(ip + OverWriteFirst(fsign[1])) * fsign[1]);
        //                 const lid_t idx2 = (!isphys[2]) ? i2 : (-(ip + OverWriteFirst(fsign[2])) * fsign[2]);

        //                 // check that we stay at the correct sport for everybody
        //                 m_assert((fstart[0] + idx0) >= (-block->gs()) && (fstart[0] + idx0) < (block->stride()), "index 0 is wrong: %d with gs = %d and stride = %d", fstart[0] + idx0, block->gs(), block->stride());
        //                 m_assert((fstart[1] + idx1) >= (-block->gs()) && (fstart[1] + idx1) < (block->stride()), "index 1 is wrong: %d with gs = %d and stride = %d", fstart[1] + idx1, block->gs(), block->stride());
        //                 m_assert((fstart[2] + idx2) >= (-block->gs()) && (fstart[2] + idx2) < (block->stride()), "index 2 is wrong: %d with gs = %d and stride = %d", fstart[2] + idx2, block->gs(), block->stride());
        //                 // check the specific isphys direction
        //                 m_assert((((fstart[0] + idx0) * isphys[0]) >= 0) && ((fstart[0] + idx0) * isphys[0]) < (block->stride() - block->gs()), "index 0 is wrong: %d with gs = %d and stride = %d", fstart[0] + idx0, block->gs(), block->stride());
        //                 m_assert((((fstart[1] + idx1) * isphys[1]) >= 0) && ((fstart[1] + idx1) * isphys[1]) < (block->stride() - block->gs()), "index 1 is wrong: %d with gs = %d and stride = %d", fstart[1] + idx1, block->gs(), block->stride());
        //                 m_assert((((fstart[2] + idx2) * isphys[2]) >= 0) && ((fstart[2] + idx2) * isphys[2]) < (block->stride() - block->gs()), "index 2 is wrong: %d with gs = %d and stride = %d", fstart[2] + idx2, block->gs(), block->stride());

        //                 // store the result
        //                 f[ip] = ldata[m_midx(idx0, idx1, idx2, 0, block)];
        //                 // get the data position
        //                 real_t data_pos[3];
        //                 m_pos_relative(data_pos, idx0, idx1, idx2, hgrid);
        //                 xf[ip] = data_pos[dir];
        //             }
        //             // get the position of the ghost
        //             real_t pos[3];
        //             m_pos_relative(pos, i0, i1, i2, hgrid);

        //             // get the ghost value
        //             ldata[m_midx(i0, i1, i2, 0, block)] = Stencil_(f, xf, pos[dir], boundary_condition);
        //         }
        //     }
        // }
    }
    //-------------------------------------------------------------------------
};

class ZeroBoundary : public Boundary<0> {
   protected:
    // real_t Stencil_(const real_t* f, const real_t x, const real_t h, const real_t normal, const real_t value) override;
    inline real_t Stencil_(m_ptr<const real_t> f, m_ptr<const real_t> xf, const real_t xp, const real_t bc_value)  override { return 0.0; };
};

template <lda_t len>
real_t NevilleInterp(real_t* c, real_t* d, const real_t* x, const real_t xp) {
    //-------------------------------------------------------------------------
    real_t value = c[0];  // ns = 0
    // m = 1 -> npoint
    for (lda_t m = 1; m < len; ++m) {
        for (lda_t i = 0; i < (len - m); ++i) {
            const real_t ho  = x[i] - xp;
            const real_t hp  = x[i + m] - xp;
            const real_t den = x[i] - x[i + m];
            const real_t w   = (c[i + 1] - d[i]) / den;
            d[i]             = hp * w;
            c[i]             = ho * w;
        }
        value += c[0];  // ns = 0
    }
    return value;
    //-------------------------------------------------------------------------
}

template <lda_t npoint>
class ExtrapBoundary : public Boundary<npoint> {
   protected:
    /**
    * @brief do not overwrite the first element if the face is oriented negatively (sign = -1)
    */
    virtual inline bool OverWriteFirst(const real_t sign) const override { return (sign > 0.5); }

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
    inline real_t Stencil_(m_ptr<const real_t> f, m_ptr<const real_t> xf, const real_t xp, const real_t bc_value)  override {
        //-------------------------------------------------------------------------
        constexpr lda_t len = npoint;
        // get the arrays
        real_t d[len];
        real_t c[len];
        real_t x[len];

        // m = 0
        for (lda_t id = 0; id < npoint; ++id) {
            d[id] = f()[id];
            c[id] = f()[id];
            x[id] = xf()[id];
        }

        // return the value
        return NevilleInterp<len>(c, d, x, xp);
        //-------------------------------------------------------------------------
    }
};

template <lda_t npoint>
class DirichletBoundary : public Boundary<npoint> {
   protected:
    /**
     * @brief Impletement the Neville algorithm to obtain an ODD polynomial around a value, i.e. a Dirichlet boundary condition
     * 
     * We refer to Numerical Recipes in C book for more details.
     * Compared to the version presented p109, we always have an ns = 0 as we always "extrapolate on the 0th side of the f array"
     * 
     * -> this one is easy, we know the value to use for the boundary
     */
    inline real_t Stencil_(m_ptr<const real_t> f, m_ptr<const real_t> xf, const real_t xp, const real_t bc_value)  override {
        m_assert(bc_value == 0.0, "to be honest I have never tested with nn-zero bc_value...");
        //-------------------------------------------------------------------------
        constexpr lda_t len = npoint + 1;
        // get the arrays
        real_t d[len];
        real_t c[len];
        real_t x[len];

        // set the bc value
        c[0] = bc_value;
        d[0] = bc_value;
        x[0] = 0.0;
        // set the rest of the vector
        for (lda_t id = 0; id < npoint; ++id) {
            // store the value
            c[id + 1] = f()[id];
            d[id + 1] = f()[id];
            x[id + 1] = xf()[id];
        }

        // return the value
        return NevilleInterp<len>(c, d, x, xp);
        //-------------------------------------------------------------------------
    }
};

template <lda_t npoint>
class NeumanBoundary : public Boundary<npoint> {
   protected:
    /**
     * @brief Impletement the Neville algorithm to obtain an EVEN polynomial around a flux, i.e. a Neuman boundary condition
     * 
     * We refer to Numerical Recipes in C book for more details.
     * Compared to the version presented p109, we always have an ns = 0 as we always "extrapolate on the 0th side of the f array"
     * 
     * Here, we have to retrieve the value of the boundary.
     * We know that with P+1 points, we can build of polynomial of order P and obtain it's derivative exactly.
     * So by inverting the FD stencil, we can match the BC and obtain the value to set
     * 
     * @param f 
     * @param xf 
     * @param xp 
     * @param bc_value is the flux value COMING IN * hgrid !!!!
     * @return real_t 
     */
    inline real_t Stencil_(m_ptr<const real_t> f, m_ptr<const real_t> xf, const real_t xp, const real_t bc_value)  override {
        //-------------------------------------------------------------------------
        m_assert(bc_value == 0.0, "to be honest I have never tested with nn-zero bc_value...");
        // get the correct
        real_t f0_value;
        if constexpr (npoint == 1) {
            f0_value = (bc_value - (1.0 * f()[0])) * (-1.0);
        } else if (npoint == 3) {
            f0_value = (bc_value - (3.0 * f()[0] - 3.0 / 2.0 * f()[1] + 1.0 / 3.0 * f()[2])) * (-6.0 / 11.0);
        } else if (npoint == 5) {
            f0_value = (bc_value - (5.0 * f()[0] - 5.0 * f()[1] + 10.0 / 3.0 * f()[2] - 5.0 / 4.0 * f()[3] + 1.0 / 5.0 * f()[4])) * (-60.0 / 137.0);
        } else {
            m_assert(false, "error, the npoint = %d is not valid", npoint);
        }
        //-------------------------------------------------------------------------
        constexpr lda_t len = npoint + 1;
        // get the arrays
        real_t d[len];
        real_t c[len];
        real_t x[len];

        // set the bc value
        c[0] = f0_value;
        d[0] = f0_value;
        x[0] = 0.0;

        // set the rest of the vector
        for (lda_t id = 0; id < npoint; ++id) {
            // store the value
            c[id + 1] = f()[id];
            d[id + 1] = f()[id];
            x[id + 1] = xf()[id];
        }

        return NevilleInterp<len>(c, d, x, xp);
        //-------------------------------------------------------------------------
    }
};

#endif  // SRC_BOUNDARY_HPP_
