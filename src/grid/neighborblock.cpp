#include "grid/neighborblock.hpp"
#include "core/types.hpp"
#include "core/macros.hpp"

inline void face_sign(const bool mask_bool, const iface_t iface, iface_t* face_dir, real_t sign[3]) {
    //-------------------------------------------------------------------------
    const real_t  mask = static_cast<real_t>(mask_bool);
    const iface_t dir  = mask_bool ? (iface / 2) : 0;
    m_assert((0 <= dir) && (dir < 3), "the dir must be between 0 and 2 and not %d", dir);
    // store the dir
    (*face_dir) += (mask_bool) ? (dir) : 0;
    // update the sign
    sign[dir] += mask * (((iface % 2) == 1) ? 1.0 : -1.0);
    // sign[(dir + 1) % 3] += 0.0;
    // sign[(dir + 2) % 3] += 0.0;
    //-------------------------------------------------------------------------
}

inline void edge_sign(const bool mask_bool, const iface_t iedge, iface_t* edge_dir, real_t sign[3]) {
    /*
    the plane convention for the sign variable convention for the sign
    2 +--------------+ 3
      |              |
      |              |
      |dir2          |
      |              |
    0 +--------------+ 1
        dir1
    */
    //-------------------------------------------------------------------------
    const real_t  mask = static_cast<real_t>(mask_bool);
    const iface_t dir  = mask_bool ? (iedge / 4) : 0;         // this is the direction of the edge
    const iface_t dir1 = static_cast<iface_t>(0 == dir);      // dir1 in the plane: dir1 = x if dir = y or z, or y if dir = x
    const iface_t dir2 = 2 - static_cast<iface_t>(2 == dir);  // dir2 in the plane: dir2 = y if dir=z, = z if dir=x or dir = y
    m_assert((0 <= dir) && (dir < 3), "the dir must be between 0 and 2 and not %d", dir);
    m_assert((0 <= dir1) && (dir1 < 3), "the dir must be between 0 and 2 and not %d", dir1);
    m_assert((0 <= dir2) && (dir2 < 3), "the dir must be between 0 and 2 and not %d", dir2);
    m_assert((dir != dir1) && (dir1 != dir2) && (dir != dir2), "the dir must be differents: %d %d %d", dir, dir1, dir2);
    // store the dir
    (*edge_dir) += (mask_bool) ? (dir) : 0;
    // update the sign
    // sign[dir] += 0.0;
    sign[dir1] += mask * (((iedge % 4) % 2) == 1 ? +1.0 : -1.0);
    sign[dir2] += mask * (((iedge % 4) / 2) == 1 ? +1.0 : -1.0);
    //-------------------------------------------------------------------------
}

inline void corner_sign(const bool mask_bool, const iface_t icorner, real_t sign[3]) {
    //-------------------------------------------------------------------------
    const real_t mask = static_cast<real_t>(mask_bool);
    sign[0] += mask * ((icorner % 2) == 1 ? +1.0 : -1.0);
    sign[1] += mask * (((icorner % 4) / 2) == 1 ? +1.0 : -1.0);
    sign[2] += mask * ((icorner / 4) == 1 ? +1.0 : -1.0);
    //-------------------------------------------------------------------------
}

/**
 * @brief Given a face, edge or a corner returns the outgoing normal (sign)
 *
 * @param ibidule for a face (`ibidule<6`), an edge (`6<= ibidule < 18`) or a corner (`18<= ibidule < 26`). can also be -1, we return 0.0 then
 * @param sign the sign of the outgoing normal
 */
void GhostGetSign(const iface_t ibidule, real_t sign[3]) {
    //-------------------------------------------------------------------------
    iface_t dir = 0;
    sign[0]     = 0.0;
    sign[1]     = 0.0;
    sign[2]     = 0.0;
    m_assert(sign[0] == 0.0 && sign[1] == 0.0 && sign[2] == 0.0, "wrong sign value: %e %e %e", sign[0], sign[1], sign[2]);

    // let's go
    const bool face_mask = (0 <= ibidule) && (ibidule < 6);
    face_sign(face_mask, ibidule, &dir, sign);
    const bool edge_mask = (6 <= ibidule) && (ibidule < 18);
    edge_sign(edge_mask, ibidule - 6, &dir, sign);
    const bool corner_mask = (18 <= ibidule) && (ibidule < 26);
    corner_sign(corner_mask, ibidule - 18, sign);

    m_assert((sign[0] == 0.0) || (sign[0] == 1.0) || (sign[0] == -1.0), "wrong sign value: %e", sign[0]);
    m_assert((sign[1] == 0.0) || (sign[1] == 1.0) || (sign[1] == -1.0), "wrong sign value: %e", sign[1]);
    m_assert((sign[2] == 0.0) || (sign[2] == 1.0) || (sign[2] == -1.0), "wrong sign value: %e", sign[2]);
    //-------------------------------------------------------------------------
};