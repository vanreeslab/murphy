#ifndef SRC_TOOLSP4EST_HPP_
#define SRC_TOOLSP4EST_HPP_

#include <limits>
#include <list>

#include "murphy.hpp"
#include "p8est.h"
#include "p8est_bits.h"
#include "p8est_ghost.h"
#include "p8est_mesh.h"

using std::numeric_limits;

template <typename T>
inline static T p4est_GetElement(sc_array_t* array, const int id) {
    return *(reinterpret_cast<T*>(sc_array_index_int(array, id)));
};

template <typename T>
inline static T* p4est_GetPointer(sc_array_t* array, const int id) {
    return reinterpret_cast<T*>(sc_array_index_int(array, id));
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

inline static int p4est_GetOwnerFromGhost(p8est_t* forest, p8est_quadrant_t* ghost) {
    //---------------------------------------------------------------------
    p4est_topidx_t tree_id = ghost->p.piggy3.which_tree;
    return p8est_quadrant_find_owner(forest, tree_id, -1, ghost);
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

inline static real_t p4est_QuadLen(const level_t level) {
    //---------------------------------------------------------------------
    const real_t val = 1.0 / (P8EST_ROOT_LEN / P8EST_QUADRANT_LEN(level));
    m_assert(val == m_quad_len(level), "euuh");
    return val;
    //---------------------------------------------------------------------
}

inline static int p4est_GetChildID(const real_t xyz[3], const level_t level) {
    //---------------------------------------------------------------------
    // mimic the behavior of p8est_quadrant_child_id (p4est_bits.c)
    const real_t len_coarse = p4est_QuadLen(level - 1);
    int          id         = 0;
    id += (fmod(xyz[0], len_coarse) == 0.0) ? 0 : 1;
    id += (fmod(xyz[1], len_coarse) == 0.0) ? 0 : 2;
    id += (fmod(xyz[2], len_coarse) == 0.0) ? 0 : 4;
    return id;
    //---------------------------------------------------------------------
}

inline static void p4est_GetNeighbor(/* p4est arguments */ p8est_t* forest, p8est_connectivity_t* connect, p8est_ghost_t* ghost, p8est_mesh_t* mesh,
                                     /* looking for */ const p4est_topidx_t tree_id, const p4est_locidx_t local_id, const iface_t ibidule,
                                     /* result */ std::list<qdrt_t*>* ngh_list, std::list<rank_t>* rank_list) {
    //---------------------------------------------------------------------
    // clear the lists
    ngh_list->clear();
    rank_list->clear();

    // get me as quad
    p8est_tree_t*     my_tree = p8est_tree_array_index(forest->trees, tree_id);
    p8est_quadrant_t* quad    = p8est_quadrant_array_index(&my_tree->quadrants, local_id);

    rank_t my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    //................................................
    // if (ibidule < 18) {
    // this is a face or an edge, we are good to go old style
    // temporary sc array used to get the ghosts
    sc_array_t* ngh_quad = sc_array_new(sizeof(p8est_quadrant_t*));  // points to the quad
    sc_array_t* ngh_qid  = sc_array_new(sizeof(int));                // give the ID of the quad or ghost
    sc_array_t* ngh_enc  = sc_array_new(sizeof(int));                // get the status

    // get the ghost
    p4est_locidx_t cum_id = my_tree->quadrants_offset + local_id;
    p8est_mesh_get_neighbors(forest, ghost, mesh, cum_id, ibidule, ngh_quad, ngh_enc, ngh_qid);

    // add them to the list
    for (iblock_t ib = 0; ib < ngh_quad->elem_count; ++ib) {
        // p8est_quadrant_t* ngh = p8est_quadrant_array_index(ngh_quad, ib);
        p8est_quadrant_t* ngh = p4est_GetElement<p8est_quadrant_t*>(ngh_quad, ib);

        // finish the check
        m_verb("adding %p to the list",ngh);
        ngh_list->push_back(ngh);

        // find out the rank, if not ghost, myself, if ghost, take the piggy number
        const int  status  = p4est_GetElement<int>(ngh_enc, ib);
        const bool isghost = (status < 0);

        if (!isghost) {
            m_verb("pushing to list: rank %d", my_rank);
            rank_list->push_back(my_rank);
        } else {
            m_verb("found block at level %d, local num = %d and tree %d", ngh->level, ngh->p.piggy3.local_num, ngh->p.piggy3.which_tree);
            const int ngh_rank = p8est_quadrant_find_owner(forest, ngh->p.piggy3.which_tree, -1, ngh);
            m_verb("pushing to list: rank %d", ngh_rank);
            rank_list->push_back(ngh_rank);
        }
    }
    
    // if it's a corner, do an extra check with the new method
    if (ibidule >= 18) {
        // we are looking for a corner
        // initialize the arrays
        sc_array_t* treeid_list    = sc_array_new(sizeof(p4est_topidx_t));
        sc_array_t* surrogate_list = sc_array_new(sizeof(p8est_quadrant_t));
        sc_array_t* exist_arr      = sc_array_new(sizeof(int));
        sc_array_t* rank_arr       = sc_array_new(sizeof(int));
        sc_array_t* quad_arr       = sc_array_new(sizeof(p8est_quadrant_t));

        // get the surrogate block
        const iface_t corner_id = ibidule - 18;
        p8est_quadrant_corner_neighbor_extra(quad, tree_id, corner_id, surrogate_list, treeid_list, NULL, connect);
        m_verb("looking for block @tree %d  and corner %d", tree_id, corner_id);
        m_verb("we found %d surrogate blocks", surrogate_list->elem_count);

        // for each surrogate block
        for (sid_t is = 0; is < surrogate_list->elem_count; is++) {
            // get the quad and the tree info

            // p8est_quadrant_t* quad_2_find = p8est_quadrant_array_index(surrogate_list, is);
            // p4est_topidx_t    tree_2_find = *(p4est_topidx_t*)sc_array_index(treeid_list, is);
            p8est_quadrant_t* quad_2_find = p4est_GetPointer<p8est_quadrant_t>(surrogate_list, is);
            p4est_topidx_t    tree_2_find = p4est_GetElement<p4est_topidx_t>(treeid_list, is);

            // reset the lists
            sc_array_reset(exist_arr);
            sc_array_reset(rank_arr);
            sc_array_reset(quad_arr);

            // actually search for it
            int is_valid = p8est_quadrant_exists(forest, ghost, tree_2_find, quad_2_find, exist_arr, rank_arr, quad_arr);
            m_assert(quad_arr->elem_count == 1, "there is %ld quad matching the needed one", quad_arr->elem_count);
            m_verb("the quadrant found is valid? %d", is_valid);
            m_verb("we have %d elements in the exist vector", exist_arr->elem_count);
            // for (int i = 0; i < exist_arr->elem_count; ++i) {
            //     m_log("value: %d", p4est_GetElement<int>(exist_arr, i));
            // }
            // m_assert(exist_arr->elem_count == 1, "there is more than 1 quad matching the needed one");
            // m_assert(*((int*)sc_array_index_int(exist_arr, 0)) == 1, "the quadrant must exist");

            // store the quadrant and the rank -> use the piggy3 to get the correct quad
            p8est_quadrant_t* quad_piggy   = p4est_GetPointer<p8est_quadrant_t>(quad_arr, 0);
            p4est_topidx_t    rank_to_push = p4est_GetElement<int>(rank_arr, 0);
            if (is_valid && rank_to_push == my_rank) {
                // we are not ghost use the piggy3 to recover the real quad
                m_verb("piggy3: tree id = %d, local_num = %d", quad_piggy->p.piggy3.which_tree, quad_piggy->p.piggy3.local_num);
                p8est_tree_t*     tree_to_push = p8est_tree_array_index(forest->trees, quad_piggy->p.piggy3.which_tree);
                p8est_quadrant_t* quad_to_push = p8est_quadrant_array_index(&tree_to_push->quadrants, quad_piggy->p.piggy3.local_num);

                if (ngh_quad->elem_count > 0) {
                    m_assert(quad_to_push == ngh_list->back(), "the quad should be the same...: %p vs %p", quad_to_push, ngh_list->back());
                } else {
                    ngh_list->push_back(quad_to_push);
                    rank_list->push_back(rank_to_push);
                }

                m_verb("pushing to list: adress: %p  and rank %d", quad_to_push, rank_to_push);

            } else if (is_valid) {
                // we are a ghost, we can push the piggy quad
                // search for the quad in the ghost then
                p8est_tree_t* tree_to_push = p8est_tree_array_index(forest->trees, quad_piggy->p.piggy3.which_tree);
                ssize_t       ghost_offset = p8est_ghost_bsearch(ghost, rank_to_push, quad_piggy->p.which_tree, quad_piggy);
                p8est_quadrant_t* ghost_to_push = p8est_quadrant_array_index(&ghost->ghosts,ghost_offset);

                if (ngh_quad->elem_count > 0) {
                    m_assert(ghost_to_push == ngh_list->back(), "the quad should be the same...: %p vs %p, piggy 3 = %d %d vs %d %d", quad_piggy, ngh_list->back(),
                             ghost_to_push->p.piggy3.which_tree, ghost_to_push->p.piggy3.local_num, ngh_list->back()->p.piggy3.which_tree, ngh_list->back()->p.piggy3.local_num);
                } else {
                    ngh_list->push_back(ghost_to_push);
                    rank_list->push_back(rank_to_push);
                }
            }
        }

        // kill the arrays
        sc_array_destroy(surrogate_list);
        sc_array_destroy(treeid_list);
        sc_array_destroy(exist_arr);
        sc_array_destroy(rank_arr);
        sc_array_destroy(quad_arr);
    }

    // destroy the arrays
    sc_array_destroy(ngh_quad);
    sc_array_destroy(ngh_qid);
    sc_array_destroy(ngh_enc);

    //---------------------------------------------------------------------
    m_assert(ngh_list->size() == rank_list->size(), "the arrays must have the same length");
}

#endif  // SRC_TOOLSP4EST_HPP_