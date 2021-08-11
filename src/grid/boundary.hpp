#ifndef SRC_GRID_BOUNDARY_HPP_
#define SRC_GRID_BOUNDARY_HPP_

#include "core/forloop.hpp"
#include "core/macros.hpp"
// #include "core/pointers.hpp"
#include "grid/cartblock.hpp"
#include "grid/subblock.hpp"

/**
 * @brief functor that implements a boundary conditions that uses some points inside the block to extrapolate a ghost value
 * 
 * @tparam the number of points used within the block to build the extrapolation
 */
template <lda_t npoint>
class Boundary {
   protected:
    /**
     * @brief Compute the a polynomial ghost value based on known field values and a boundary condition
     * 
     * @param x_f the array of the position of known field values (do not contain the BC)
     * @param f the know field values (do not contain the BC)
     * @param xb the position of the BC
     * @param fb the BC
     * @param xp the position of the wanted value
     * @return real_t the value of the polynomial extrapolated/interpolated to that point
     * 
     * if @ref OverWriteFirst_() is true:
     * ```
     *          GHOST   |                       BLOCK                      |    GHOST
     * -----------------|--------------------------------------------------|---------------
     *      normal = -1 |                                                  | normal = +1
     *                  |                                                  |
     *           x (<0) |                                                  |   x (>0) 
     *     <----------->|   xf[0]   xf[1]                 xf[1]  xf[0]     |<-------->
     * ---o------o------o------o------o---   ...       -----o------o-------o------o------
     *                       f[0]   f[1]     ...         f[1]    f[0]
     * ```
     * 
     * if @ref OverWriteFirst_() is false:
     * ```
     *          GHOST   |                       BLOCK                      |    GHOST
     * -----------------|--------------------------------------------------|---------------
     *      normal = -1 |                                                  | normal = +1
     *                  |                                                  |
     *           x (<0) |                                                  |   x (>0) 
     *     <----------->|   xf[1]   xf[2]                 xf[2]  xf[1]     |<-------->
     * ---o------o------o------o------o---   ...       -----o------o-------o------o------
     *                 f[0]   f[1]   f[2]     ...         f[2]    f[1]    f[0]
     * ```
     */
    // virtual inline real_t Stencil_(const real_t* const  f, const real_t* const  xf, const real_t xp, const real_t bc_value) const = 0;
    virtual inline real_t Stencil_(const real_t* const xf, const real_t* const f, const real_t xb, const real_t fb, const real_t xp) const = 0;

    /**
     * @brief returns true if the border element of the block is overwritten
     * 
     * by default we overwrite the last element (fsign > 0 ) but not the first one (fsign < 0)
     */
    virtual inline bool OverWriteFirst_(const real_t fsign) const { return (fsign > (0.5)); };

   public:
    /**
     * @brief compute the boundary condition on a block, in the area of a ghost block
     * 
     * @param iface the index of the face on which we apply it (X- = 0, X+ = 1, Y- = 2, Y+ = 3, Z- = 4, Z- = 5)
     * @param first_last the index of the first/last point in the block
     * @param boundary_condition the BC value
     * @param block the original block (needed for hgrid and block start/end)
     * @param gblock the ghost block where to apply it
     * @param data the data
     */
    virtual void operator()(const sid_t iface, const bidx_t first_last[3], const real_t hgrid[3], const real_t boundary_condition, const SubBlock* const gblock, const data_ptr data) const {
        m_assert(iface >= 0 && iface < 6, "iface = %d is wrong", iface);
        m_assert((npoint + 1) <= M_N, "the size of the box is not big enough to take the needed samples");
        // m_assert(block->stride() == gblock->stride(), "the two strides must be the same");
        // m_assert(block->gs() == gblock->gs(), "the two gs must be the same");
        //-------------------------------------------------------------------------
        const real_t face_sign[6][3] = {{-1.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, -1.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, -1.0}, {0.0, 0.0, 1.0}};

        // get the face direction and other informations
        const iface_t dir       = iface / 2;
        const bool    isphys[3] = {dir == 0, dir == 1, dir == 2};
        const real_t* fsign     = face_sign[iface];

        // move the data around the starting point of the face to fill
        bidx_t  b_stride = gblock->stride();
        real_t* sdata    = data.Write(0);
        m_assert(gblock->stride() == data.stride() && gblock->gs() == data.gs(), "The layout and the memory must be equal");

        // m_log("sdata from data = %ld", (sdata - data()) / sizeof(real_t));
        // m_log("the boundary will be from %d %d %d to %d %d %d", gblock->start(0), gblock->start(1), gblock->start(2), gblock->end(0), gblock->end(1), gblock->end(2));

        auto op = [=, &sdata](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
            // get the local information
            real_t* ldata = sdata + m_idx(i0, i1, i2, 0, b_stride);

            // get the distance between me and the face
            // the face_start is given as M_N, so if I am positive signed, I need add 1 to reach the last element of the block
            const bidx_t dis0 = i0 - first_last[0];  //+ 2 * (fsign[0] > 0.5);
            const bidx_t dis1 = i1 - first_last[1];  //+ 2 * (fsign[1] > 0.5);
            const bidx_t dis2 = i2 - first_last[2];  //+ 2 * (fsign[2] > 0.5);

            // m_log("local is %d %d %d -> dir = %d -> hgrid = %e ", i0, i1, i2,dir,hgrid[dir]);
            // m_log("distance is %d %d %d", dis0, dis1, dis2);
            // m_log("ghost stride = %d", b_stride);

            // we need interpolations points in the face direction
            real_t f[npoint];
            real_t xf[npoint];
            for (bidx_t ip = 0; ip < npoint; ip++) {
                // if the direction is not physics, just consider the current ID,
                // if the direction is physics, takes the first, second and third point INSIDE the block (0 = first inside)
                // these are local increments on the source which is already in position fstart
                const bidx_t idx0 = (!isphys[0]) ? 0 : (-(ip + OverWriteFirst_(fsign[0])) * fsign[0] - dis0);
                const bidx_t idx1 = (!isphys[1]) ? 0 : (-(ip + OverWriteFirst_(fsign[1])) * fsign[1] - dis1);
                const bidx_t idx2 = (!isphys[2]) ? 0 : (-(ip + OverWriteFirst_(fsign[2])) * fsign[2] - dis2);

                // check that we stay at the correct sport for everybody
                m_assert((i0 + idx0) >= (-gblock->gs()) && (i0 + idx0) < (b_stride), "index 0 is wrong: %d with gs = %d and stride = %d", i0 + idx0, gblock->gs(), b_stride);
                m_assert((i1 + idx1) >= (-gblock->gs()) && (i1 + idx1) < (b_stride), "index 1 is wrong: %d with gs = %d and stride = %d", i1 + idx1, gblock->gs(), b_stride);
                m_assert((i2 + idx2) >= (-gblock->gs()) && (i2 + idx2) < (b_stride), "index 2 is wrong: %d with gs = %d and stride = %d", i2 + idx2, gblock->gs(), b_stride);
                // check the specific isphys direction
                m_assert((((i0 + idx0) * isphys[0]) >= 0) && ((i0 + idx0) * isphys[0]) < (b_stride - gblock->gs()), "index 0 is wrong: %d with gs = %d and stride = %d", i0 + idx0, gblock->gs(), gblock->stride());
                m_assert((((i1 + idx1) * isphys[1]) >= 0) && ((i1 + idx1) * isphys[1]) < (b_stride - gblock->gs()), "index 1 is wrong: %d with gs = %d and stride = %d", i1 + idx1, gblock->gs(), gblock->stride());
                m_assert((((i2 + idx2) * isphys[2]) >= 0) && ((i2 + idx2) * isphys[2]) < (b_stride - gblock->gs()), "index 2 is wrong: %d with gs = %d and stride = %d", i2 + idx2, gblock->gs(), gblock->stride());

                // store the result
                f[ip] = ldata[m_idx(idx0, idx1, idx2, 0, b_stride)];

                // get the data position, relative to me, i.e. using the idx indexes
                // const real_t data_pos[3] = {(idx0)*hgrid[0],
                //                             (idx1)*hgrid[1],
                //                             (idx2)*hgrid[2]};
                // xf[ip] = data_pos[dir];

                // xf[ip] = idx0 * (dir == 0) * hgrid[0] + idx1 * (dir == 1) * hgrid[1] + idx2 * (dir == 2) * hgrid[2];
                xf[ip] = idx0 * (dir == 0) + idx1 * (dir == 1) + idx2 * (dir == 2);
            }
            // get the ghost value
            // ldata[0] = Stencil_(f, xf, 0.0, boundary_condition);
            const real_t xb = -dis0 * (dir == 0) - dis1 * (dir == 1) - dis2 * (dir == 2);
            // const real_t xb = -dis0 * (dir == 0) * hgrid[0] - dis1 * (dir == 1) * hgrid[1] - dis2 * (dir == 2) * hgrid[2];
            ldata[0]        = Stencil_(xf, f, xb, boundary_condition, 0.0);
        };

        // m_log("doing boundary from %d %d %d to %d %d %d", gblock->start(0), gblock->start(1), gblock->start(2), gblock->end(0), gblock->end(1), gblock->end(2));
        // get the starting and ending indexes
        const bidx_t start[3] = {gblock->start(0), gblock->start(1), gblock->start(2)};
        const bidx_t end[3]   = {gblock->end(0), gblock->end(1), gblock->end(2)};
        for_loop(&op, start, end);
    }
    //-------------------------------------------------------------------------
};

class ZeroBoundary : public Boundary<0> {
   protected:
    inline bool   OverWriteFirst_(const real_t fsign) const override { return true; }
    inline real_t Stencil_(const real_t* const xf, const real_t* const f, const real_t xb, const real_t fb, const real_t xp) const override { return 0.0; };
};

template <lda_t len>
real_t NevilleInterp(real_t* c, real_t* d, const real_t* x, const real_t xp, short_t ns) {
    //-------------------------------------------------------------------------
    real_t value = c[ns--];

    // m = 1 -> npoint
    for (short_t m = 1; m < len; ++m) {
        for (short_t i = 0; i < (len - m); ++i) {
            const real_t ho  = x[i] - xp;
            const real_t hp  = x[i + m] - xp;
            const real_t den = ho - hp;  //x[i] - x[i + m];
            const real_t w   = (c[i + 1] - d[i]);
            d[i]             = (hp / den) * w;
            c[i]             = (ho / den) * w;
        }
        // value += c[0];  // ns = 0
        value += (2 * (ns + 1) < (len - m)) ? c[ns + 1] : d[ns--];
    }
    return value;
    //-------------------------------------------------------------------------
}

/**
 * @brief Extrapolating boundary condition
 * 
 * @tparam npoint polynomial order = npoint - 1
 */
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
    // inline real_t Stencil_(const real_t* const  f, const real_t* const  xf, const real_t xp, const real_t bc_value, const rea) const override {
    inline real_t Stencil_(const real_t* const xf, const real_t* const f, const real_t xb, const real_t fb, const real_t xp) const override {
        //-------------------------------------------------------------------------
        constexpr lda_t len = npoint;
        // get the arrays
        real_t d[len];
        real_t c[len];
        real_t x[len];

        // m = 0
        for (lda_t id = 0; id < len; ++id) {
            c[id] = f[id];
            d[id] = f[id];
            x[id] = xf[id];
        }

        // return the value
        return NevilleInterp<len>(c, d, x, xp, 0);
        //-------------------------------------------------------------------------
    }
};

/**
 * @brief Dirichlet boundary conditions
 * 
 * @tparam npoint polynomial order = npoint
 */
template <lda_t npoint>
class DirichletBoundary : public Boundary<npoint> {
   protected:
    /**
    * @brief do overwrite the first element
    */
    virtual inline bool OverWriteFirst_(const real_t fsign) const override { return true; }

    /**
     * @brief Impletement the Neville algorithm to obtain an ODD polynomial around a value, i.e. a Dirichlet boundary condition
     * 
     * We refer to Numerical Recipes in C book for more details.
     * Compared to the version presented p109, we always have an ns = 0 as we always "extrapolate on the 0th side of the f array"
     * 
     * -> this one is easy, we know the value to use for the boundary
     */
    // inline real_t Stencil_(const real_t* const  f, const real_t* const  xf, const real_t xp, const real_t bc_value) const override {
    inline real_t Stencil_(const real_t* const xf, const real_t* const f, const real_t xb, const real_t fb, const real_t xp) const override {
        m_assert(fb == 0.0, "to be honest I have never tested with nn-zero bc_value...");
        //-------------------------------------------------------------------------
        constexpr lda_t len = npoint + 1;
        // get the arrays
        real_t d[len];
        real_t c[len];
        real_t x[len];

        // set the bc value
        c[0] = fb;
        d[0] = fb;
        x[0] = xb;
        // set the rest of the vector
        for (lda_t id = 0; id < npoint; ++id) {
            // store the value
            c[id + 1] = f[id];
            d[id + 1] = f[id];
            x[id + 1] = xf[id];
        }

        // return the value
        return NevilleInterp<len>(c, d, x, xp, 0);
        //-------------------------------------------------------------------------
    }
};

/**
 * @brief Neuman boundary conditions
 * 
 * @tparam npoint polynomial order = npoint
 */
template <lda_t npoint>
class NeumanBoundary : public Boundary<npoint> {
   protected:
    /**
    * @brief do overwrite the first element
    */
    virtual inline bool OverWriteFirst_(const real_t fsign) const override { return true; }

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
    // inline real_t Stencil_(const real_t* const  f, const real_t* const  xf, const real_t xp, const real_t bc_value) const override {
    inline real_t Stencil_(const real_t* const xf, const real_t* const f, const real_t xb, const real_t fb, const real_t xp) const override {
        m_assert(fb == 0.0, "to be honest I have never tested with nn-zero bc_value...");
        //-------------------------------------------------------------------------

        // get the correct
        real_t f0_value;
        if constexpr (npoint == 1) {
            f0_value = (fb - (1.0 * f[0])) * (-1.0);
        } else if (npoint == 3) {
            f0_value = (fb - (3.0 * f[0] - 3.0 / 2.0 * f[1] + 1.0 / 3.0 * f[2])) * (-6.0 / 11.0);
        } else if (npoint == 5) {
            f0_value = (fb - (5.0 * f[0] - 5.0 * f[1] + 10.0 / 3.0 * f[2] - 5.0 / 4.0 * f[3] + 1.0 / 5.0 * f[4])) * (-60.0 / 137.0);
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
        x[0] = xb;

        // set the rest of the vector
        for (lda_t id = 0; id < npoint; ++id) {
            // store the value
            c[id + 1] = f[id];
            d[id + 1] = f[id];
            x[id + 1] = xf[id];
        }

        return NevilleInterp<len>(c, d, x, xp, 0);
        //-------------------------------------------------------------------------
    }
};

#endif  // SRC_BOUNDARY_HPP_
