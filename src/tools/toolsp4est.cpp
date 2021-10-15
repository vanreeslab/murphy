
#include "toolsp4est.hpp"

/**
 * @brief return lists with: the quadrants, the cummulative id and the rank to access the neighor
 * 
 * @param local_id the local (non-cummulative!!) id of the quadrant
 * 
 */
void p4est_GetNeighbor(/* p4est arguments */ p8est_t* forest, p8est_connectivity_t* connect, p8est_ghost_t* ghost, p8est_mesh_t* mesh,
                       /* looking for */ const p4est_topidx_t tree_id, const p4est_locidx_t local_id, const iface_t ibidule,
                       /* result */ std::list<qdrt_t*>* ngh_list, std::list<iblock_t>* id_list, std::list<rank_t>* rank_list) {
    //--------------------------------------------------------------------------
    // clear the lists
    ngh_list->clear();
    id_list->clear();
    rank_list->clear();

    // get me as quad
    p8est_tree_t*     my_tree = p8est_tree_array_index(forest->trees, tree_id);
    p8est_quadrant_t* quad    = p8est_quadrant_array_index(&my_tree->quadrants, local_id);

    rank_t my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    //................................................
    // we always try first the good old way and correct afterwards if needed
    // temporary sc array used to get the ghosts
    sc_array_t* ngh_quad = sc_array_new(sizeof(p8est_quadrant_t*));  // points to the quad
    sc_array_t* ngh_qid  = sc_array_new(sizeof(int));                // give the ID of the quad or ghost
    sc_array_t* ngh_enc  = sc_array_new(sizeof(int));                // get the status

    // get the ghost
    p4est_locidx_t cum_id = my_tree->quadrants_offset + local_id;
    p8est_mesh_get_neighbors(forest, ghost, mesh, cum_id, ibidule, ngh_quad, ngh_enc, ngh_qid);

    // add them to the list
    for (iblock_t ib = 0; ib < ngh_quad->elem_count; ++ib) {
        m_assert(ngh_qid->elem_count == ngh_quad->elem_count, "the counters = %ld and %ld must be =", ngh_qid->elem_count, ngh_quad->elem_count);
        m_assert(ngh_enc->elem_count == ngh_quad->elem_count, "the counters = %ld and %ld must be =", ngh_enc->elem_count, ngh_quad->elem_count);
        m_verb("ngh_qid = %d, ngh_enc = %d, ngh_quad = %d", ngh_qid->elem_count, ngh_quad->elem_count, ngh_enc->elem_count);
        // p8est_quadrant_t* ngh = p8est_quadrant_array_index(ngh_quad, ib);

        p8est_quadrant_t* ngh = p4est_GetElement<p8est_quadrant_t*>(ngh_quad, ib);
        ngh_list->push_back(ngh);

        // find out the rank, if not ghost, myself, if ghost, take the piggy number
        const int  status  = p4est_GetElement<int>(ngh_enc, ib);
        const bool isghost = (status < 0);
        m_verb("block %d searching for ghost #%d: status = %d", cum_id, ibidule, status);

        if (!isghost) {
            m_verb("pushing to list: rank %d", my_rank);
            rank_list->push_back(my_rank);
            // if it's a local block, we trust the id
            int ngh_block_id = p4est_GetElement<int>(ngh_qid, ib);
            m_assert(ngh_block_id < std::numeric_limits<iblock_t>::max(), "the id %d must be smaller than the limit %d ", ngh_block_id, std::numeric_limits<iblock_t>::max());
            id_list->push_back(ngh_block_id);
        } else {
            m_verb("found block at level %d, local num = %d and tree %d", ngh->level, ngh->p.piggy3.local_num, ngh->p.piggy3.which_tree);
            const int ngh_rank = p8est_quadrant_find_owner(forest, ngh->p.piggy3.which_tree, -1, ngh);
            m_verb("pushing to list: rank %d", ngh_rank);
            rank_list->push_back(ngh_rank);

            // if it's a ghost we use the piggy3 id instead (they don't match, no clue why)
            int ngh_block_id = ngh->p.piggy3.local_num;  // p4est_GetElement<int>(ngh_qid, ib);
            m_assert(ngh_block_id < std::numeric_limits<iblock_t>::max(), "the id %d must be smaller than the limit %d ", ngh_block_id, std::numeric_limits<iblock_t>::max());
            id_list->push_back(ngh_block_id);
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

        // for each surrogate block
        for (sid_t is = 0; is < surrogate_list->elem_count; is++) {
            // get the quad and the tree info
            p8est_quadrant_t* quad_2_find = p4est_GetPointer<p8est_quadrant_t>(surrogate_list, is);
            p4est_topidx_t    tree_2_find = p4est_GetElement<p4est_topidx_t>(treeid_list, is);

            // reset the lists
            sc_array_reset(exist_arr);
            sc_array_reset(rank_arr);
            sc_array_reset(quad_arr);

            // actually search for it
            int is_valid = p8est_quadrant_exists(forest, ghost, tree_2_find, quad_2_find, exist_arr, rank_arr, quad_arr);
            m_assert(quad_arr->elem_count == 1, "there is %ld quad matching the needed one", quad_arr->elem_count);
            m_verb("the quadrant found is valid? %d", is_valid);  // i don't understand this value....
            m_verb("we have %d elements in the exist vector", exist_arr->elem_count);

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

                    // get the cummulative ID
                    int ngh_block_id = quad_piggy->p.piggy3.local_num;
                    m_assert(ngh_block_id < std::numeric_limits<iblock_t>::max(), "the id %d must be smaller than the limit %d ", ngh_block_id, std::numeric_limits<iblock_t>::max());
                    id_list->push_back(ngh_block_id);
                }

                m_verb("pushing to list: adress: %p  and rank %d", quad_to_push, rank_to_push);
            } else if (is_valid) {
                // we are a ghost, we can push the piggy quad
                // search for the quad in the ghost then
                p8est_tree_t*     tree_to_push  = p8est_tree_array_index(forest->trees, quad_piggy->p.piggy3.which_tree);
                bidx_t            ghost_offset  = p8est_ghost_bsearch(ghost, rank_to_push, quad_piggy->p.which_tree, quad_piggy);
                p8est_quadrant_t* ghost_to_push = p8est_quadrant_array_index(&ghost->ghosts, ghost_offset);

                if (ngh_quad->elem_count > 0) {
                    m_assert(ghost_to_push == ngh_list->back(), "the quad should be the same...: %p vs %p, piggy 3 = %d %d vs %d %d", quad_piggy, ngh_list->back(),
                             ghost_to_push->p.piggy3.which_tree, ghost_to_push->p.piggy3.local_num, ngh_list->back()->p.piggy3.which_tree, ngh_list->back()->p.piggy3.local_num);
                } else {
                    ngh_list->push_back(ghost_to_push);
                    rank_list->push_back(rank_to_push);
                    // get the cummulative ID
                    int ngh_block_id = ghost_to_push->p.piggy3.local_num;
                    m_assert(ngh_block_id < std::numeric_limits<iblock_t>::max(), "the id %d must be smaller than the limit %d ", ngh_block_id, std::numeric_limits<iblock_t>::max());
                    id_list->push_back(ngh_block_id);
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

    //--------------------------------------------------------------------------
    m_assert(ngh_list->size() == rank_list->size(), "the arrays must have the same length");
}
