#ifndef SRC_TOOLSP4EST_HPP_
#define SRC_TOOLSP4EST_HPP_

#include <limits>
#include <list>

#include "core/macros.hpp"
#include "core/types.hpp"
#include "p8est.h"
#include "p8est_bits.h"
#include "p8est_ghost.h"
#include "p8est_mesh.h"

using std::numeric_limits;

struct p4est_Essentials {
    const bool*    is_periodic;
    p8est_t*       forest;
    p8est_mesh_t*  mesh;
    p8est_ghost_t* ghost;
};

template <typename T>
inline static T p4est_GetElement(sc_array_t* array, const int id) {
    return *(reinterpret_cast<T*>(sc_array_index_int(array, id)));
};

template <typename T>
inline static T* p4est_GetPointer(sc_array_t* array, const int id) {
    return reinterpret_cast<T*>(sc_array_index_int(array, id));
};

inline static char p4est_MaxLocalLevel(const p8est_t* forest) {
    //--------------------------------------------------------------------------
    char l_max_level = 0;
    for (p4est_topidx_t it = forest->first_local_tree; it <= forest->last_local_tree; it++) {
        // get the current tree
        p8est_tree_t* ctree = p8est_tree_array_index(forest->trees, it);
        // get the max delta level over the current tree
        l_max_level = m_max(ctree->maxlevel, l_max_level);  // max level is given by the tree
    }
    return l_max_level;
    //--------------------------------------------------------------------------
};

inline static int p4est_GetOwnerFromGhost(p8est_t* forest, p8est_quadrant_t* ghost) {
    //--------------------------------------------------------------------------
    p4est_topidx_t tree_id = ghost->p.piggy3.which_tree;
    return p8est_quadrant_find_owner(forest, tree_id, -1, ghost);
    //--------------------------------------------------------------------------
};

inline static p8est_quadrant_t* p4est_GetQuadFromMirror(const p8est_t* forest, const p8est_quadrant_t* mirror) {
    //--------------------------------------------------------------------------
    p8est_tree_t*  tree    = p8est_tree_array_index(forest->trees, mirror->p.piggy3.which_tree);
    p4est_locidx_t quad_id = mirror->p.piggy3.local_num - tree->quadrants_offset;
    return p8est_quadrant_array_index(&tree->quadrants, quad_id);
    //--------------------------------------------------------------------------
};

inline static p4est_locidx_t p4est_NumQuadOnLevel(const p8est_mesh_t* mesh, const char level) {
    m_assert(level >= 0, "the level = %d must be >=0", level);
    m_assert(level < P8EST_MAXLEVEL, "the level = %d must be <= %d", level, P8EST_MAXLEVEL);
    //--------------------------------------------------------------------------
    size_t num = mesh->quad_level[level].elem_count;
    m_assert(num < numeric_limits<p4est_locidx_t>::max(), "the number of element is too big to be local");
    return num;
    //--------------------------------------------------------------------------
};
inline static p4est_locidx_t p4est_GetQuadIdOnLevel(const p8est_mesh_t* mesh, const char level, const p4est_locidx_t quad_id) {
    m_assert(quad_id < numeric_limits<int>::max(), "quad id is too big");
    //--------------------------------------------------------------------------
    sc_array_t quad_id_array = mesh->quad_level[level];
    return p4est_GetElement<p4est_locidx_t>(&quad_id_array, (int)quad_id);
    //--------------------------------------------------------------------------
};

inline static real_t p4est_QuadLen(const level_t level) {
    m_assert(level >= 0, "the level = %d must be >=0", level);
    //--------------------------------------------------------------------------
    // const real_t val = 1.0 / (real_t)(P8EST_ROOT_LEN / P8EST_QUADRANT_LEN(level));
    const real_t val = 1.0 / ((real_t)(1 << (level)));
    m_assert(val > 0.0, "the length = %e must be >0", val);
    return val;
    //--------------------------------------------------------------------------
};

inline static int p4est_GetChildID(const real_t xyz[3], const level_t level) {
    m_assert(level > 0, "the level = %d must be > 0", level);
    //--------------------------------------------------------------------------
    // mimic the behavior of p8est_quadrant_child_id (p4est_bits.c)
    const real_t len_coarse = p4est_QuadLen(level - 1);

    int id = 0;
    id += (fmod(xyz[0], len_coarse) == 0.0) ? 0 : 1;
    id += (fmod(xyz[1], len_coarse) == 0.0) ? 0 : 2;
    id += (fmod(xyz[2], len_coarse) == 0.0) ? 0 : 4;
    return id;
    //--------------------------------------------------------------------------
};

/**
 * @brief gets the Block of type B stored as a user_data of the p4est quadrant
 */
template <typename B>
B* p4est_GetBlock(const qdrt_t* quad) {
    //-------------------------------------------------------------------------------
    B** p4est_usr_data = static_cast<B**>(quad->p.user_data);
    B*  block          = p4est_usr_data[0];
    m_assert(block != nullptr, "the block address cannot be null");
    return block;
    //-------------------------------------------------------------------------------
}

/**
 * @brief sets the Block of type B stored as a user_data of the p4est quadrant
 */
template <typename B>
static inline void p4est_SetBlock(qdrt_t* quad, B* block) {
    //-------------------------------------------------------------------------------
    m_assert(block != nullptr, "the block address cannot be null, we have an issue here");
    B** p4est_usr_data = static_cast<B**>(quad->p.user_data);
    p4est_usr_data[0]  = block;
    //-------------------------------------------------------------------------------
}

void p4est_GetNeighbor(/* p4est arguments */ p8est_t* forest, p8est_connectivity_t* connect, p8est_ghost_t* ghost, p8est_mesh_t* mesh,
                       /* looking for */ const p4est_topidx_t tree_id, const p4est_locidx_t local_id, const iface_t ibidule,
                       /* result */ std::list<qdrt_t*>* ngh_list, std::list<iblock_t>* id_list, std::list<rank_t>* rank_list);

#endif  // SRC_TOOLSP4EST_HPP_
