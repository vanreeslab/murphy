#ifndef SRC_TOOLSP4EST_HPP_
#define SRC_TOOLSP4EST_HPP_

#include "p8est.h"
#include "p8est_mesh.h"
#include <limits>

using std::numeric_limits;

template <typename T>
inline static T p4est_GetElement(sc_array_t* array, const int id) {
    return *(reinterpret_cast<T*>(sc_array_index_int(array, id)));
};

inline static int8_t p4est_MaxLocalLevel(const p8est_t* forest) {
    //---------------------------------------------------------------------
    int8_t l_max_level = 0;
    for (p4est_topidx_t it = forest->first_local_tree; it <= forest->last_local_tree; it++) {
        // get the current tree
        p8est_tree_t* ctree = p8est_tree_array_index(forest->trees, it);
        // get the max delta level over the current tree
        l_max_level = m_max(ctree->maxlevel, l_max_level);  // max level is given by the tree
    }
    return l_max_level;
    //---------------------------------------------------------------------
};

inline static p8est_quadrant_t* p4est_GetQuadFromMirror(const p8est_t* forest, p8est_quadrant_t* mirror) {
    //---------------------------------------------------------------------
    p8est_tree_t*  tree    = p8est_tree_array_index(forest->trees, mirror->p.piggy3.which_tree);
    p4est_locidx_t quad_id = mirror->p.piggy3.local_num - tree->quadrants_offset;
    return p8est_quadrant_array_index(&tree->quadrants, quad_id);
    //---------------------------------------------------------------------
};

inline static p4est_locidx_t p4est_NumQuadOnLevel(const p8est_mesh_t* mesh, const int8_t level) {
    //---------------------------------------------------------------------
    size_t num = mesh->quad_level[level].elem_count;
    m_assert(num < numeric_limits<p4est_locidx_t>::max(), "the number of element is too big to be local");
    return num;
    //---------------------------------------------------------------------
};
inline static p4est_locidx_t p4est_GetQuadIdOnLevel(const p8est_mesh_t* mesh, const int8_t level, const p4est_locidx_t quad_id) {
    m_assert(quad_id < numeric_limits<int>::max(), "quad id is too big");
    //---------------------------------------------------------------------
    sc_array_t quad_id_array = mesh->quad_level[level];
    return p4est_GetElement<p4est_locidx_t>(&quad_id_array, (int)quad_id);
    //---------------------------------------------------------------------
};

#endif  // SRC_TOOLSP4EST_HPP_