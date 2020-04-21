#include "boundary.hpp"

static const real_t face_sign[6][3]  = {{-1.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, -1.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, -1.0}, {0.0, 0.0, 1.0}};
static const lid_t  face_start[6][3] = {{0, 0, 0}, {M_N - 1, 0, 0}, {0, 0, 0}, {0, M_N - 1, 0}, {0, 0, 0}, {0, 0, M_N - 1}};


void Boundary::operator()(const real_t boundary_condition, const real_t hgrid[3], PhysBlock* block, real_p data) {
    const lid_t*  fstart = face_start[block->iface()];
    const real_t* fsign  = face_sign[block->iface()];
    // get the face direction and a boolean array with it
    const sid_t dir       = block->iface() / 2;
    const bool  isphys[3] = {dir == 0, dir == 1, dir == 2};
    // get the offset in memory, i.e. the position of the first point (0,0,0)
    real_t offset[3];
    m_pos_relative(offset,0,0,0,hgrid);

    // shift the data to the correct spot
    real_p ldata = data + m_idx(fstart[0], fstart[1], fstart[2]);

    // let's goooo
    for (lid_t i2 = block->start(2); i2 < block->end(2); i2++) {
        for (lid_t i1 = block->start(1); i1 < block->end(1); i1++) {
            for (lid_t i0 = block->start(0); i0 < block->end(0); i0++) {
                // we need three interpolations points in the face direction
                real_t f[3];

                // if the direction is not physics, just consider the current ID,
                // if the direction is physics, takes the first, second and third point INSIDE the block (0 = first inside)
                f[0] = ldata[m_idx((!isphys[0]) * i0,
                                   (!isphys[1]) * i1,
                                   (!isphys[2]) * i2)];
                f[1] = ldata[m_idx((!isphys[0]) * i0 - 1 * fsign[0] * isphys[0],
                                   (!isphys[1]) * i1 - 1 * fsign[1] * isphys[1],
                                   (!isphys[2]) * i2 - 1 * fsign[2] * isphys[2])];
                f[2] = ldata[m_idx((!isphys[0]) * i0 - 2 * fsign[0] * isphys[0],
                                   (!isphys[1]) * i1 - 2 * fsign[1] * isphys[1],
                                   (!isphys[2]) * i2 - 2 * fsign[2] * isphys[2])];
                // get the position
                const real_t x = (i0 * isphys[0] + i1 * isphys[1] + i2 * isphys[2]) * hgrid[dir];
                // get the ghost value
                ldata[m_idx(i0, i1, i2)] = Stencil_(f, x, hgrid[dir], offset[dir], fsign[dir], boundary_condition);
            }
        }
    }
}

real_t EvenBoundary_4::Stencil_(const real_t* f, const real_t x, const real_t h, const real_t offset, const real_t normal, const real_t flux) {
    // remove the flux slope form the value
    // if x is positive, we need go in the negative and vice-versa
    real_t cor_f[3] = {f[0] - normal * flux * (0.0 + offset) * h,
                       f[1] - normal * flux * (1.0 + offset) * h,
                       f[2] - normal * flux * (2.0 + offset) * h};
    // compute the stencil coefficients
    const real_t a = 1.0 / (12.0 * pow(h, 4)) * f[0] - 1.0 / (8.0 * pow(h, 4)) * f[1] + 1.0 / (24.0 * pow(h, 4)) * f[2];
    const real_t b = -17.0 / (24.0 * pow(h, 2)) * f[0] + 13.0 / (16.0 * pow(h, 2)) * f[1] - 5.0 / (48.0 * pow(h, 2)) * f[2];
    const real_t c = (75.0 / 64.0) * f[0] - (25.0 / 128.0) * f[1] + (3.0 / 128.0) * f[2];
    // return the polynomial + the flux slope
    return a * pow(x, 4) + b * pow(x, 2) + c + flux * x;
}
