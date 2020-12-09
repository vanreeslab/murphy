#include "ghost.hpp"

#include <algorithm>
#include <limits>

#include "boundary.hpp"
#include "defs.hpp"
#include "omp.h"
#include "toolsp4est.hpp"
#include "wavelet.hpp"

#define M_NNEIGHBOR 26

/**
 * @brief the localization of the interface
 */
static lid_t face_start[6][3] = {{0, 0, 0}, {M_N, 0, 0}, {0, 0, 0}, {0, M_N, 0}, {0, 0, 0}, {0, 0, M_N}};

/**
 * @brief Construct a new Ghost object 
 * 
 * see @ref Ghost::Ghost(ForestGrid* grid, const level_t min_level, const level_t max_level, Wavelet* interp) for details
 * 
 * @param grid the ForestGrid to use, must have been initiated using ForestGrid::SetupP4estGhostMesh() 
 * @param interp the Wavelet to use, will drive the number of ghost points to consider
 */
Ghost::Ghost(ForestGrid* grid, const Wavelet* interp, Prof* profiler) : Ghost(grid, -1, P8EST_MAXLEVEL + 1, interp, profiler) {
    //-------------------------------------------------------------------------
    // we called the function Ghost::Ghost(ForestGrid* grid, const level_t min_level, const level_t max_level, Wavelet* interp)
    //-------------------------------------------------------------------------
}

/**
 * @brief Construct a new Ghost, allocate the ghost lists and initiates the communications
 * 
 * Once created, the ghost is fixed for a given grid. if the grid changes, a new Ghost objects has to be created.
 * 
 * @note: The Wavelet has to be given beforehands because the it drives the number of actual GP to consider.
 * While for memory alignement, the number of ghost points is given by M_GS, the wavelet does not require that many ghost points to
 * be computed. To reduce the memory cost, only the needed ghost points will be computed.
 * 
 * @param grid the ForestGrid to use, must have been initiated using ForestGrid::SetupP4estGhostMesh() 
 * @param min_level the minimum level on which the GP are initiated
 * @param max_level the maximum level on which the GP are initiated
 * @param interp the Wavelet to use, will drive the number of ghost points to consider
 */
Ghost::Ghost(ForestGrid* grid, const level_t min_level, const level_t max_level, const Wavelet* interp, Prof* profiler) : interp_(interp) {
    m_begin;
    m_assert(grid->is_mesh_valid(), "the mesh needs to be valid before entering here");
    //-------------------------------------------------------------------------
    // get the important pointers
    grid_   = grid;
    interp_ = interp;
    prof_   = profiler;
    // store the level information
    min_level_ = m_max(min_level, 0);
    max_level_ = m_min(max_level, P8EST_QMAXLEVEL);

    //................................................
    // get how many active quads should be considered
    // n_active_quad_     = 0;
    // for (level_t il = min_level_; il <= max_level_; il++) {
    //     n_active_quad_ += p4est_NumQuadOnLevel(mesh, il);
    // }
    // m_verb("I will ghost %d local active quads", n_active_quad_);

    // store the number of ghosts needed
    m_assert(interp_->nghost_front() <= M_GS, "The memory for the ghost points is too small: M_GS = %d vs nghost = %d", M_GS, interp_->nghost_front());
    m_assert(interp_->nghost_back() <= M_GS, "The memory for the ghost points is too small: M_GS = %d vs nghost = %d", M_GS, interp_->nghost_back());

    //................................................
    // initialize the communications and the ghost's lists
    m_profStart(prof_, "ghost_init");
    InitComm_();
    InitList_();
    m_profStop(prof_, "ghost_init");

    //-------------------------------------------------------------------------
    m_log("ghost for refinement: %d %d", interp_->nrefine_front(), interp_->nrefine_back());
    m_log("ghost for coarsening: %d %d", interp_->ncoarsen_front(), interp_->ncoarsen_back());
    m_log("ghost for criterion: %d %d + shift: %d %d", interp_->ncriterion_front(), interp_->ncriterion_back(), interp_->shift_front(), interp_->shift_back());
    m_log("ghost initialized with %s, nghost = %d %d, coarse nghost = %d %d", interp_->Identity().c_str(), interp_->nghost_front(), interp_->nghost_back(), interp_->CoarseNGhostFront(), interp_->CoarseNGhostBack());
    m_end;
}

/**
 * @brief Destroy the Ghost and free the allocated memory
 * 
 */
Ghost::~Ghost() {
    m_begin;
    m_assert(grid_->is_mesh_valid(), "the mesh needs to be valid before entering here");
    //-------------------------------------------------------------------------
    FreeList_();
    FreeComm_();
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief initialize the list of ghost blocks: every block knows where to find the ghosts afterwards
 * 
 */
void Ghost::InitList_() {
    m_begin;
    m_verb("Ghost lists initialization started...");
    //-------------------------------------------------------------------------
    // get stupid MPI info
    rank_t         mpi_size = grid_->mpisize();
    MPI_Comm       mpi_comm = grid_->mpicomm();
    p8est_t*       forest   = grid_->p4est_forest();
    p8est_ghost_t* ghost    = grid_->p4est_ghost();
    p8est_mesh_t* mesh = grid_->p4est_mesh();

    // sanity checks
    m_assert(mpi_comm == MPI_COMM_WORLD, "the comm should be a comm world");

    //................................................
    // allocate the array to link the local_id to the mirror displacement and init it
    MPI_Info info;
    MPI_Info_create(&info);
    MPI_Info_set(info, "no_locks", "true");
    MPI_Aint win_mem_size = mesh->local_num_quadrants * sizeof(MPI_Aint);
    // check the size
    m_assert(win_mem_size >= 0, "the memory size should be >=0");
    // allocate the array and the window with the offsets
    local2disp_ = reinterpret_cast<MPI_Aint*>(m_calloc(win_mem_size));
    MPI_Win_create(local2disp_, win_mem_size, sizeof(MPI_Aint), info, mpi_comm, &local2disp_window_);
    // MPI_Win_allocate(win_mem_size, sizeof(MPI_Aint), info, mpi_comm, &local2disp_, &local2disp_window_);
    MPI_Info_free(&info);
    m_verb("allocating %ld bytes in the window for %d active quad", win_mem_size, mesh->local_num_quadrants);
    // make sure everything is done
    // MPI_Win_fence(0, local2disp_window_);
    MPI_Barrier(MPI_COMM_WORLD);

    //................................................
    // compute the number of admissible local mirrors and store their reference in the array
    iblock_t active_mirror_count = 0;
    for (iblock_t im = 0; im < ghost->mirrors.elem_count; im++) {
        qdrt_t* mirror       = p8est_quadrant_array_index(&ghost->mirrors, im);
        level_t mirror_level = p4est_GetQuadFromMirror(forest, mirror)->level;
        // update the counters if the mirror is admissible
        // we need to have called the function p4est_balance()!!
        if ((min_level_ - 1) <= mirror_level && mirror_level <= (max_level_ + 1)) {
            // store the displacement in the local2disp_ array
            iblock_t local_id     = mirror->p.piggy3.local_num;
            local2disp_[local_id] = (active_mirror_count++) * m_blockmemsize(1);
        }
    }
    // make sure everybody did it
    // MPI_Win_fence(0, local2disp_window_);
    MPI_Barrier(MPI_COMM_WORLD);

    //................................................
    // post the exposure epoch and start the access one for local2mirrors
    m_assert(mirror_origin_group_ != MPI_GROUP_NULL, "call the InitComm function first!");
    m_assert(mirror_target_group_ != MPI_GROUP_NULL, "call the InitComm function first!");

    // start the exposure epochs if any
    MPI_Win_post(mirror_origin_group_, 0, local2disp_window_);
    // we need to not start the access epochs if we are empty (it fails on the beast otherwise)
    if (mirror_target_group_ != MPI_GROUP_EMPTY) {
        MPI_Win_start(mirror_target_group_, 0, local2disp_window_);
    }
    // check that we will NOT go for GhostInitLists
    m_assert(!(active_mirror_count != 0 && mirror_target_group_ == MPI_GROUP_EMPTY && mpi_size > 1), "if we have some active quadrant, we need to start the get epochs: active_mirror_count = %d", active_mirror_count);

    // init the list on every active block that matches the level requirements
    // const Wavelet*
    const ForestGrid* mygrid = grid_;
    for (level_t il = min_level_; il <= max_level_; il++) {
        // DoOpMeshLevel(this, &Ghost::InitList4Block_, grid_, il);
        DoOpMeshLevel(nullptr, &GridBlock::GhostInitLists, mygrid, il, mygrid, interp_, local2disp_window_);
    }

    // complete the epoch and wait for the exposure one
    if (mirror_target_group_ != MPI_GROUP_EMPTY) {
        MPI_Win_complete(local2disp_window_);
    }
    // stop the exposure epochs if any
    MPI_Win_wait(local2disp_window_);

    //................................................
    MPI_Win_free(&local2disp_window_);
    m_free(local2disp_);
    local2disp_window_ = MPI_WIN_NULL;
    local2disp_        = nullptr;

    //-------------------------------------------------------------------------
    m_verb("Ghost lists initialization is done");
    m_end;
}

void Ghost::FreeList_() {
    m_begin;
    //-------------------------------------------------------------------------
    for (level_t il = min_level_; il <= max_level_; il++) {
        // DoOpMeshLevel(this, &Ghost::FreeList4Block_, grid_, il);
        DoOpMeshLevel(nullptr, &GridBlock::GhostFreeLists, grid_, il);
    }
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief Initialize the communication pattern: init the windows for RMA access
 * 
 * We initiate the window mirrors_window_ and the associated groups: mirror_origin_group_ & mirror_target_group_
 * 
 */
void Ghost::InitComm_() {
    m_begin;
    // m_assert(n_active_quad_ >= 0, "the number of active quads must be computed beforehand");
    //-------------------------------------------------------------------------
    // get stupid information
    int            mpi_size = grid_->mpisize();
    MPI_Comm       mpi_comm = grid_->mpicomm();
    p8est_t*       forest   = grid_->p4est_forest();
    p8est_ghost_t* ghost    = grid_->p4est_ghost();

    //................................................
    // compute the number of admissible local mirrors and store their reference in the array
    n_mirror_to_send_ = 0;
    for (iblock_t im = 0; im < ghost->mirrors.elem_count; im++) {
        qdrt_t* mirror       = p8est_quadrant_array_index(&ghost->mirrors, im);
        level_t mirror_level = p4est_GetQuadFromMirror(forest, mirror)->level;
        // update the counters if the mirror is admissible (i.e. it satisfies the requirements)
        // we need to have called the function p4est_balance()!!
        if ((min_level_ - 1) <= mirror_level && mirror_level <= (max_level_ + 1)) {
            n_mirror_to_send_++;
        }
    }
    m_verb("I have %d mirrors to send", n_mirror_to_send_);
    // initialize the Window by allocating the memory space needed for the exchange
    MPI_Info info;
    MPI_Info_create(&info);
    MPI_Info_set(info, "no_locks", "true");
    MPI_Aint win_mem_size = n_mirror_to_send_ * m_blockmemsize(1) * sizeof(real_t);
    mirrors_              = reinterpret_cast<real_t*>(m_calloc(win_mem_size));
    MPI_Win_create(mirrors_, win_mem_size, sizeof(real_t), info, mpi_comm, &mirrors_window_);
    MPI_Info_free(&info);
    // because of alignement issues (https://github.com/open-mpi/ompi/issues/7955), cannot use this one
    // MPI_Win_allocate(win_mem_size, sizeof(real_t), info, mpi_comm, &mirrors_, &mirrors_window_);
    m_assert(m_isaligned(mirrors_), "the mirror temp array is not aligned");
    m_assert(mirrors_window_ != MPI_WIN_NULL, "the MPI window created is null, which is not a good news");

    //................................................
    // get the list of ranks that will generate a call to access my mirrors
    rank_t n_in_group = 0;
    // over allocate the array to its max size and fill it only partially... not great
    rank_t* group_ranks = reinterpret_cast<rank_t*>(m_calloc(mpi_size * sizeof(rank_t)));
    for (rank_t ir = 0; ir < mpi_size; ir++) {
        // for every mirror send to the current
        iblock_t send_first = ghost->mirror_proc_offsets[ir];
        iblock_t send_last  = ghost->mirror_proc_offsets[ir + 1];
        for (iblock_t bid = send_first; bid < send_last; bid++) {
            // access the element and ask for its level
            qdrt_t* mirror       = p8est_quadrant_array_index(&ghost->mirrors, ghost->mirror_proc_mirrors[bid]);
            level_t mirror_level = p4est_GetQuadFromMirror(forest, mirror)->level;
            // if the mirror is admissible, register the rank and break
            if ((min_level_ - 1) <= mirror_level && mirror_level <= (max_level_ + 1)) {
                group_ranks[n_in_group++] = ir;
                break;
            }
        }
    }
    // get the RMA mirror group ready - the group that will need my info
    MPI_Group win_group;
    MPI_Win_get_group(mirrors_window_, &win_group);
    MPI_Group_incl(win_group, n_in_group, group_ranks, &mirror_origin_group_);

    //................................................
    // get the list of ranks that will received a call from me to access their mirrors
    n_in_group = 0;
    for (rank_t ir = 0; ir < mpi_size; ir++) {
        // for the current rank, determine if I have ghosts to him
        iblock_t recv_first = ghost->proc_offsets[ir];
        iblock_t recv_last  = ghost->proc_offsets[ir + 1];
        for (iblock_t bid = recv_first; bid < recv_last; bid++) {
            // access the element and ask for its level
            p8est_quadrant_t* ghostquad   = p8est_quadrant_array_index(&ghost->ghosts, bid);
            level_t           ghost_level = ghostquad->level;
            // if the ghost block is admissible, register the rank and break
            if ((min_level_ - 1) <= ghost_level && ghost_level <= (max_level_ + 1)) {
                group_ranks[n_in_group++] = ir;
                break;
            }
        }
    }
    // add the cpus that will get a call from me
    MPI_Group_incl(win_group, n_in_group, group_ranks, &mirror_target_group_);
    MPI_Group_free(&win_group);

    //................................................
    // free the allocated memory
    m_free(group_ranks);

    // assert that everybody has created the windows correctly
    // MPI_Win_fence(0, mirrors_window_);
    MPI_Barrier(MPI_COMM_WORLD);
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief Free the communication window mirrors_window_ and the groups mirror_origin_group_ & mirror_target_group_
 * 
 */
void Ghost::FreeComm_() {
    m_begin;
    //-------------------------------------------------------------------------
    // free the group
    MPI_Group_free(&mirror_origin_group_);
    MPI_Group_free(&mirror_target_group_);
    // free the window
    m_free(mirrors_);
    MPI_Win_free(&mirrors_window_);
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief Post the RMA call: get the values of siblings and coarse neighbors to my local memory
 * 
 * - push the current info to the windows so that other ranks can access it
 * - init the epochs on the windows
 * - call Ghost::GetGhost4Block_Post()
 * 
 * @param field 
 * @param ida 
 */
void Ghost::PullGhost_Post(const Field* field, const lda_t ida) {
    m_begin;
    m_assert(ida >= 0, "the ida must be >=0!");
    m_assert(grid_->is_mesh_valid(), "the mesh needs to be valid before entering here");
    //-------------------------------------------------------------------------
    // store the current dimension
    ida_ = ida;

    //................................................
    // fill the Window memory with the Mirror information
    m_profStart(prof_, "ghost computation");
    LoopOnMirrorBlock_(&Ghost::PushToWindow4Block, field);
    m_profStop(prof_, "ghost computation");

    //................................................
    // post the exposure epoch for my own mirrors: I am a target warning that origin group will RMA me
    MPI_Win_post(mirror_origin_group_, 0, mirrors_window_);
    // start the access epoch, to get info from neighbors: I am an origin warning that I will RMA the target group
    if (mirror_target_group_ != MPI_GROUP_EMPTY) {
        MPI_Win_start(mirror_target_group_, 0, mirrors_window_);
    }

    //................................................
    // start what can be done = sibling and parents local copy + physical BC + myself copy
    m_profStart(prof_, "ghost computation");
    for (level_t il = min_level_; il <= max_level_; il++) {
        // DoOpMeshLevel(this, &Ghost::GetGhost4Block_Post, grid_, il, field);
        DoOpMeshLevel(nullptr, &GridBlock::GhostGet_Post, grid_, il, field, ida, interp_, mirrors_window_);
    }
    m_profStop(prof_, "ghost computation");
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief Wait for the RMA calls to finish and update the coarser neighbor's ghosting points
 * 
 * @param field the field on which to operate
 * @param ida the dimemsion inside the field
 */
void Ghost::PullGhost_Wait(const Field* field, const lda_t ida) {
    m_begin;
    m_assert(ida >= 0, "the ida must be >=0!");
    m_assert(ida_ == ida, "the ongoing dimension (%d) must be over first", ida_);
    m_assert(grid_->is_mesh_valid(), "the mesh needs to be valid before entering here");
    //-------------------------------------------------------------------------

    //................................................
    // finish the access epochs for the exposure epoch to be over
    m_profStart(prof_, "ghost wait");
    if (mirror_target_group_ != MPI_GROUP_EMPTY) {
        MPI_Win_complete(mirrors_window_);
    }
    MPI_Win_wait(mirrors_window_);
    m_profStop(prof_, "ghost wait");

    // we now have all the information needed to compute the ghost points in coarser blocks
    m_profStart(prof_, "ghost computation");
    for (level_t il = min_level_; il <= max_level_; il++) {
        // DoOpMeshLevel(this, &Ghost::GetGhost4Block_Wait, grid_, il, field);
        DoOpMeshLevel(nullptr, &GridBlock::GhostGet_Wait, grid_, il, field, ida, interp_);
    }
    m_profStop(prof_, "ghost computation");

    //................................................
    // post exposure and access epochs for to put the values to my neighbors
    MPI_Win_post(mirror_origin_group_, 0, mirrors_window_);
    if (mirror_target_group_ != MPI_GROUP_EMPTY) {
        MPI_Win_start(mirror_target_group_, 0, mirrors_window_);
    }

    //................................................
    // start what can be done = sibling and parents copy
    m_profStart(prof_, "ghost computation");
    for (level_t il = min_level_; il <= max_level_; il++) {
        // DoOpMeshLevel(this, &Ghost::PutGhost4Block_Post, grid_, il, field);
        DoOpMeshLevel(nullptr, &GridBlock::GhostPut_Post, grid_, il, field, ida, interp_, mirrors_window_);
    }
    m_profStop(prof_, "ghost computation");

    m_profStart(prof_, "ghost wait");
    // finish the access epochs for the exposure epoch to be over
    if (mirror_target_group_ != MPI_GROUP_EMPTY) {
        MPI_Win_complete(mirrors_window_);
    }
    MPI_Win_wait(mirrors_window_);
    m_profStop(prof_, "ghost wait");

    m_profStart(prof_, "ghost computation");
    // we copy back the missing info
    LoopOnMirrorBlock_(&Ghost::PullFromWindow4Block, field);

    // we now have all the information needed, we finish with a physbc
    for (level_t il = min_level_; il <= max_level_; il++) {
        // DoOpMeshLevel(this, &Ghost::PutGhost4Block_Wait, grid_, il, field);
        DoOpMeshLevel(nullptr, &GridBlock::GhostPut_Wait, grid_, il, field, ida, interp_);
    }
    m_profStop(prof_, "ghost computation");

    //-------------------------------------------------------------------------
    m_end;
}

// /**
//  * @brief setup the ghosting lists, for the considered block
//  *
//  * @warning this function is called within a omp parellel region.
//  * The p4est arrays are NOT thread-safe, while the lists have been assumed to be NOT threadsafe as well.
//  *
//  * @note the p4est array allocation is thread-safe, while the tracking of the memory allocated is NOT.
//  *
//  * @param qid the quadrant ID
//  * @param block the associated block
//  */
// void Ghost::InitList4Block_(const qid_t* qid, GridBlock* block) const {
//     // //-------------------------------------------------------------------------
//     // // get the current lists
//     // ListGBLocal*  bsibling   = block->local_sibling();
//     // ListGBLocal*  bparent    = block->local_parent();
//     // ListGBLocal*  bparent_rv = block->local_parent_reverse();
//     // ListGBMirror* gchildren  = block->ghost_children();
//     // ListGBMirror* gsibling   = block->ghost_sibling();
//     // ListGBMirror* gparent    = block->ghost_parent();
//     // ListGBMirror* gparent_rv = block->ghost_parent_reverse();
//     // listGBPhysic* phys       = block->phys();

//     // //................................................
//     // // allocate the ghost pointer
//     // block->AllocateCoarsePtr(interp_->CoarseMemSize());

//     // //................................................
//     // // temporary sc array used to get the ghosts
//     // sc_array_t* ngh_quad;  // points to the quad
//     // sc_array_t* ngh_qid;   // give the ID of the quad or ghost
//     // sc_array_t* ngh_enc;   // get the status

//     // //#pragma omp critical
//     // {
//     //     ngh_quad = sc_array_new(sizeof(qdrt_t*));
//     //     ngh_enc  = sc_array_new(sizeof(int));
//     //     ngh_qid  = sc_array_new(sizeof(int));
//     // }

//     // p8est_t*       forest = grid_->p4est_forest();
//     // p8est_mesh_t*  mesh   = grid_->p4est_mesh();
//     // p8est_ghost_t* ghost  = grid_->p4est_ghost();

//     // //................................................
//     // // get the number of ghost and the min/max of a block
//     // lid_t  nghost_front[3], nghost_back[3];
//     // lid_t  block_min[3], block_max[3];
//     // real_t block_len[3];
//     // real_t coarse_hgrid[3];
//     // for (lda_t id = 0; id < 3; id++) {
//     //     // set the number of ghost to compute
//     //     nghost_front[id] = interp_->nghost_front();
//     //     nghost_back[id]  = interp_->nghost_back();
//     //     block_min[id]    = -nghost_front[id];
//     //     block_max[id]    = M_N + nghost_back[id];
//     //     block_len[id]    = p4est_QuadLen(block->level());
//     //     coarse_hgrid[id] = p4est_QuadLen(block->level()) / M_HN;
//     // }

//     // for (iface_t ibidule = 0; ibidule < M_NNEIGHBOR; ibidule++) {
//     //     //................................................
//     //     //#pragma omp critical
//     //     {
//     //         // get the neighboring quadrant
//     //         sc_array_reset(ngh_quad);
//     //         sc_array_reset(ngh_enc);
//     //         sc_array_reset(ngh_qid);
//     //         // m_log("I have a grid with %d blocks and %d ghosts. Looking for neighbor %d of block %d = %d, %d",mesh->local_num_quadrants,ghost->ghosts.elem_count,ibidule,qid->cid,qid->qid,qid->tid);
//     //         p8est_mesh_get_neighbors(forest, ghost, mesh, qid->cid, ibidule, ngh_quad, ngh_enc, ngh_qid);
//     //     }
//     //     const iblock_t nghosts = ngh_enc->elem_count;

//     //     //................................................
//     //     // no ghosts? then is a physical BC
//     //     if (nghosts == 0) {
//     //         // we only apply the physics to entire faces
//     //         if (ibidule < 6) {
//     //             PhysBlock* pb = new PhysBlock(ibidule, block, nghost_front, nghost_back);
//     //             //#pragma omp critical
//     //             phys->push_back(pb);
//     //             m_verb("I found a physical boundary ghost!\n");
//     //         }
//     //         // else, the edges and corners will be filled through the face
//     //     }

//     //     //................................................
//     //     // this is a real block or a ghost
//     //     for (iblock_t nid = 0; nid < nghosts; nid++) {
//     //         const int  status  = *(ngh_enc->array + nid * sizeof(int));
//     //         const bool isghost = (status < 0);
//     //         qdrt_t*    nghq    = p4est_GetElement<qdrt_t*>(ngh_quad, nid);

//     //         // get the sign, i.e. the normal to the face, the edge of the corner we consider
//     //         real_t sign[3];
//     //         GhostGetSign(ibidule, sign);

//     //         //................................................
//     //         // get the position of the neighbor, as seen by me!!! may be different than the actual one if there is a periodic bc
//     //         real_t ngh_pos[3];
//     //         if (!isghost) {
//     //             // cannot use the p8est function because the which_tree is not assigned, so we retrieve the position through the block
//     //             GridBlock* ngh_block = *(reinterpret_cast<GridBlock**>(nghq->p.user_data));
//     //             ngh_pos[0]           = ngh_block->xyz(0);
//     //             ngh_pos[1]           = ngh_block->xyz(1);
//     //             ngh_pos[2]           = ngh_block->xyz(2);
//     //         } else {
//     //             p8est_qcoord_to_vertex(grid_->p4est_connect(), nghq->p.piggy3.which_tree, nghq->x, nghq->y, nghq->z, ngh_pos);
//     //         }
//     //         // fix the shift in coordinates needed IF the domain is periodic
//     //         for (lda_t id = 0; id < 3; id++) {
//     //             // if we are periodic, we overwrite the position in the direction of the normal !!ONLY!!
//     //             // since it is my neighbor in this normal direction, I am 100% sure that it's origin corresponds to the end of my block
//     //             const real_t to_replace = sign[id] * sign[id] * grid_->domain_periodic(id);  // is (+-1)^2 = +1 if we need to replace it, 0.0 otherwize
//     //             // get the expected position
//     //             const real_t expected_pos = block->xyz(id) + (sign[id] > 0.5) * p4est_QuadLen(block->level()) - (sign[id] < -0.5) * p4est_QuadLen(nghq->level);
//     //             // we override the position if a replacement is needed only
//     //             ngh_pos[id] = to_replace * expected_pos + (1.0 - to_replace) * ngh_pos[id];
//     //         }
//     //         // get the hgrid
//     //         const real_t ngh_len[3]   = {p4est_QuadLen(nghq->level), p4est_QuadLen(nghq->level), p4est_QuadLen(nghq->level)};
//     //         const real_t ngh_hgrid[3] = {p4est_QuadLen(nghq->level) / M_N, p4est_QuadLen(nghq->level) / M_N, p4est_QuadLen(nghq->level) / M_N};

//     //         //................................................
//     //         // create the new block and push back
//     //         if (!isghost) {
//     //             // associate the corresponding neighboring block
//     //             GridBlock* ngh_block = *(reinterpret_cast<GridBlock**>(nghq->p.user_data));

//     //             // register the gb in a list
//     //             if (nghq->level == block->level()) {
//     //                 // sibling: source = neighbor GridBlock, target = me
//     //                 GBLocal* gb = new GBLocal(/* source */ ngh_block->level(), ngh_pos, ngh_hgrid, ngh_len,
//     //                                           /* traget */ block->level(), block->xyz(), block->hgrid(), block_min, block_max, block->gs(), block->stride(), -1);
//     //                 gb->data_src(ngh_block);
//     //                 //#pragma omp critical
//     //                 bsibling->push_back(gb);
//     //             } else if (nghq->level < block->level()) {
//     //                 // TODO: change this to the coarse block to avoid any cast afterwards...
//     //                 // parent: source = neighbor, target = me
//     //                 GBLocal* gb = new GBLocal(/* source */ ngh_block->level(), ngh_pos, ngh_hgrid, ngh_len,
//     //                                           /* target */ block->level(), block->xyz(), block->hgrid(), block_min, block_max, block->gs(), block->stride(), -1);
//     //                 gb->data_src(ngh_block);
//     //                 //#pragma omp critical
//     //                 bparent->push_back(gb);
//     //                 // the children: the source = the coarse myself, target = my neighbor
//     //                 GBLocal* invert_gb = new GBLocal(/* source */ block->level() - 1, block->xyz(), coarse_hgrid, block_len,
//     //                                                  /* target */ ngh_block->level(), ngh_pos, ngh_hgrid, block_min, block_max, block->gs(), block->stride(), -1);
//     //                 invert_gb->data_src(ngh_block);
//     //                 //#pragma omp critical
//     //                 bparent_rv->push_back(invert_gb);
//     //             } else {
//     //                 m_assert((nghq->level - block->level()) == 1, "The delta level is not correct: %d - %d", nghq->level, block->level());
//     //             }
//     //         }
//     //         //................................................
//     //         else {
//     //             // get the local number in the remote rank and the remote rank
//     //             rank_t ngh_local_id = nghq->p.piggy3.local_num;
//     //             rank_t ngh_rank     = p4est_GetOwnerFromGhost(forest, nghq);
//     //             m_assert(ngh_rank > -1, "p4est unable to recover the rank... baaaad news");

//     //             // register the ghost block in a list
//     //             //................................................
//     //             if (nghq->level == block->level()) {
//     //                 // create the new mirror block
//     //                 // GBMirror* gb = new GBMirror(block->level(), block->xyz(), block->hgrid(), nghq->level, ngh_pos, ngh_block->hgrid(), block_len, block_min, block_max, block->gs(), block->stride(), ngh_rank);
//     //                 // sibling: source = neighbor GridBlock, target = me
//     //                 GBMirror* gb = new GBMirror(/* source */ nghq->level, ngh_pos, ngh_hgrid, ngh_len,
//     //                                             /* target */ block->level(), block->xyz(), block->hgrid(), block_min, block_max, block->gs(), block->stride(), ngh_rank);
//     //                 // ask the displacement (will be available later, when completing the call
//     //                 MPI_Get(gb->data_src_ptr(), 1, MPI_AINT, ngh_rank, ngh_local_id, 1, MPI_AINT, local2disp_window_);
//     //                 //#pragma omp critical
//     //                 gsibling->push_back(gb);
//     //             }
//     //             //................................................
//     //             else if (nghq->level < block->level()) {
//     //                 // TODO: change this to the coarse block to avoid any cast afterwards...
//     //                 // parent: source = neighbor, target = me
//     //                 GBMirror* gb = new GBMirror(/* source */ nghq->level, ngh_pos, ngh_hgrid, ngh_len,
//     //                                             /* target */ block->level(), block->xyz(), block->hgrid(), block_min, block_max, block->gs(), block->stride(), ngh_rank);
//     //                 // ask the displacement (will be available later, when completing the call
//     //                 MPI_Get(gb->data_src_ptr(), 1, MPI_AINT, ngh_rank, ngh_local_id, 1, MPI_AINT, local2disp_window_);
//     //                 //#pragma omp critical
//     //                 gparent->push_back(gb);
//     //                 // I compute my own contribution to my neighbor ghost points
//     //                 GBMirror* invert_gb = new GBMirror(/* source */ block->level() - 1, block->xyz(), coarse_hgrid, block_len,
//     //                                                    /* target */ nghq->level, ngh_pos, ngh_hgrid, block_min, block_max, block->gs(), block->stride(), ngh_rank);
//     //                 MPI_Get(invert_gb->data_src_ptr(), 1, MPI_AINT, ngh_rank, ngh_local_id, 1, MPI_AINT, local2disp_window_);
//     //                 //#pragma omp critical
//     //                 gparent_rv->push_back(invert_gb);
//     //             }
//     //             //................................................
//     //             else if (nghq->level > block->level()) {
//     //                 const real_t ngh_hgrid_coarse[3] = {ngh_len[0] / M_HN, ngh_len[1] / M_HN, ngh_len[2] / M_HN};
//     //                 // children: source = coarse version of my neighbor, target = myself
//     //                 GBMirror* gb = new GBMirror(/* source */ nghq->level - 1, ngh_pos, ngh_hgrid_coarse, ngh_len,
//     //                                             /* target */ block->level(), block->xyz(), block->hgrid(), block_min, block_max, block->gs(), block->stride(), ngh_rank);
//     //                 //#pragma omp critical
//     //                 gchildren->push_back(gb);

//     //                 m_verb("me: lvl = %d, pos = %f %f %f, length = %f", block->level(), block->xyz(0), block->xyz(1), block->xyz(2), block_len[0]);
//     //                 m_verb("ngh: lvl = %d, pos = %f %f %f, length = %f", nghq->level, ngh_pos[0], ngh_pos[1], ngh_pos[2], ngh_len[0]);
//     //                 m_verb("we have a children ghost: from %d %d %d to %d %d %d", gb->start(0), gb->start(1), gb->start(2), gb->end(0), gb->end(1), gb->end(2));
//     //                 m_verb("\n");
//     //             }
//     //             //................................................
//     //             else {
//     //                 m_assert(false, "The delta level is not correct: %d - %d", nghq->level, block->level());
//     //             }
//     //         }
//     //     }
//     // }
//     // //#pragma omp critical
//     // {
//     //     sc_array_destroy(ngh_quad);
//     //     sc_array_destroy(ngh_enc);
//     //     sc_array_destroy(ngh_qid);
//     // }
//     // //-------------------------------------------------------------------------
// }

// void Ghost::FreeList4Block_(const qid_t* qid, GridBlock* block) const {
//     //-------------------------------------------------------------------------
//     ListGBLocal*  bsibling   = block->local_sibling();
//     ListGBLocal*  bparent    = block->local_parent();
//     ListGBLocal*  bparent_rv = block->local_parent_reverse();
//     ListGBMirror* gchildren  = block->ghost_children();
//     ListGBMirror* gsibling   = block->ghost_sibling();
//     ListGBMirror* gparent    = block->ghost_parent();
//     ListGBMirror* gparent_rv = block->ghost_parent_reverse();
//     listGBPhysic* phys       = block->phys();

//     // destroy the associated blocks
//     auto remove_block = [](auto block) { delete (block); };
//     std::for_each(bsibling->begin(), bsibling->end(), remove_block);
//     std::for_each(bparent->begin(), bparent->end(), remove_block);
//     std::for_each(bparent_rv->begin(), bparent_rv->end(), remove_block);
//     std::for_each(gchildren->begin(), gchildren->end(), remove_block);
//     std::for_each(gparent->begin(), gparent->end(), remove_block);
//     std::for_each(gparent_rv->begin(), gparent_rv->end(), remove_block);
//     std::for_each(phys->begin(), phys->end(), remove_block);

//     // clear the lists
//     bsibling->clear();
//     bparent->clear();
//     bparent_rv->clear();
//     gchildren->clear();
//     gsibling->clear();
//     gparent->clear();
//     gparent_rv->clear();
//     phys->clear();
//     //-------------------------------------------------------------------------
// }

/**
 * @brief Push the local mirrors to the window's associated memory
 * 
 * @param qid the quarant id considered
 * @param block the grid block considered
 * @param fid the field ID
 */
void Ghost::PushToWindow4Block(const qid_t* qid, GridBlock* block, const Field* fid) const {
    m_assert(ida_ >= 0, "the current working dimension has to be correct");
    m_assert(ida_ < fid->lda(), "the current working dimension has to be correct");
    //-------------------------------------------------------------------------
    // recover the mirro spot using the mirror id
    real_p   mirror = mirrors_ + qid->mid * m_blockmemsize(1);
    data_ptr data   = block->pointer(fid, ida_);
    m_assume_aligned(mirror);
    m_assume_aligned(data);
    memcpy(mirror, data, m_blockmemsize(1) * sizeof(real_t));
    //-------------------------------------------------------------------------
}

/**
 * @brief Pull the local information from the window to the local block, only for the needed ghost points
 * 
 * @param qid the quarant ID considered
 * @param block the gridblock considered
 * @param fid the field id
 */
void Ghost::PullFromWindow4Block(const qid_t* qid, GridBlock* block, const Field* fid) const {
    m_assert(ida_ >= 0, "the current working dimension has to be correct");
    m_assert(ida_ < fid->lda(), "the current working dimension has to be correct");
    //-------------------------------------------------------------------------
    real_p   mirror = mirrors_ + qid->mid * m_blockmemsize(1) + m_zeroidx(0, block);
    data_ptr data   = block->data(fid, ida_);
    m_assume_aligned(mirror);
    m_assume_aligned(data);

    for (auto gblock : (*block->ghost_children())) {
        const lid_t start[3] = {gblock->start(0), gblock->start(1), gblock->start(2)};
        const lid_t end[3]   = {gblock->end(0), gblock->end(1), gblock->end(2)};

        real_t* data_src = mirror + m_midx(start[0], start[1], start[2], 0, block);
        real_t* data_trg = data + m_midx(start[0], start[1], start[2], 0, block);

        // copy the value = sendrecv to myself to the correct spot
        MPI_Status   status;
        MPI_Datatype dtype;
        ToMPIDatatype(start, end, block->stride(), 1, &dtype);
        MPI_Sendrecv(data_src, 1, dtype, 0, 0, data_trg, 1, dtype, 0, 0, MPI_COMM_SELF, &status);
        MPI_Type_free(&dtype);
    }
    //-------------------------------------------------------------------------
}

// /**
//  * @brief get the memory from the mirrors and myself (+bc!) to my ghost points and to the coarse myself
//  * 
//  * This is the step number 1/4 of ghost computation
//  * 
//  * @param qid the quarant ID considered
//  * @param block the gridblock considered
//  * @param fid the field id
//  */
// void Ghost::GetGhost4Block_Post(const qid_t* qid, GridBlock* block, const Field* fid) const {
//     //-------------------------------------------------------------------------
//     // get the working array given the thread
//     real_p     tmp       = block->coarse_ptr();
//     const bool do_coarse = (block->local_parent()->size() + block->ghost_parent()->size()) > 0;

//     //................................................
//     // start to obtain the missing info with my siblings
//     Compute4Block_GetRma2Myself_(block->ghost_sibling(), fid, block->data(fid, ida_));
//     Compute4Block_Copy2Myself_(block->local_sibling(), fid, block->data(fid, ida_));

//     //................................................
//     // if I need contributions from my parents
//     if (do_coarse) {
//         // reset the coarse memory
//         memset(tmp, 0, interp_->CoarseMemSize());

//         // get the missing part from my ghost neighbors
//         Compute4Block_GetRma2Coarse_(block->ghost_sibling(), fid, tmp);
//         Compute4Block_GetRma2Coarse_(block->ghost_parent(), fid, tmp);

//         // get the missing part from my local neighbors
//         Compute4Block_Copy2Coarse_(block->local_sibling(), fid, tmp);
//         Compute4Block_Copy2Coarse_(block->local_parent(), fid, tmp);
//     }
//     //-------------------------------------------------------------------------
// }

// /**
//  * @brief once the communication is over, compute the refined ghost points for the missing areas
//  * 
//  * This is the step number 2/4 of ghost computation
//  * 
//  * @param qid the quarant ID considered
//  * @param block the gridblock considered
//  * @param fid the field id
//  */
// void Ghost::GetGhost4Block_Wait(const qid_t* qid, GridBlock* block, const Field* fid) const {
//     //-------------------------------------------------------------------------
//     real_p     tmp       = block->coarse_ptr();
//     const bool do_coarse = (block->local_parent()->size() + block->ghost_parent()->size()) > 0;
//     if (do_coarse) {
//         // now that everything has arrived, I can compute myself and the phys boundaries
//         Compute4Block_Myself2Coarse_(qid, block, fid, tmp);
//         // refine the parents from the up-to-date coarse block
//         Compute4Block_Refine_(block->local_parent(), tmp, block->data(fid, ida_));
//         Compute4Block_Refine_(block->ghost_parent(), tmp, block->data(fid, ida_));
//     }
//     //-------------------------------------------------------------------------
// }

// /**
//  * @brief locally compute the ghost points of my coarser neighbors and post the copy/put calls
//  * 
//  * This is the step number 3/4 of ghost computation
//  * 
//  * @param qid the quarant ID considered
//  * @param block the gridblock considered
//  * @param fid the field id
//  */
// void Ghost::PutGhost4Block_Post(const qid_t* qid, GridBlock* block, const Field* fid) const {
//     //-------------------------------------------------------------------------
//     const bool do_coarse = (block->local_parent()->size() + block->ghost_parent()->size()) > 0;
//     if (do_coarse) {
//         // reset the tmp to use for the put operations
//         real_p tmp = block->coarse_ptr();
//         memset(tmp, 0, interp_->CoarseMemSize());
//         // apply the physics to the best of my knowledge
//         Compute4Block_Phys2Myself_(qid, block, fid);
//         // get the coarse representation
//         Compute4Block_Coarsen2Coarse_(block->data(fid, ida_), tmp);
//         // start the put commands
//         Compute4Block_PutRma2Parent_(block->ghost_parent_reverse(), tmp);
//         Compute4Block_Copy2Parent_(block->local_parent_reverse(), tmp, fid);
//     }
//     //-------------------------------------------------------------------------
// }

// /**
//  * @brief Finally apply the physical BC once everything is done
//  * 
//  * This is the step number 4/4 of ghost computation
//  * 
//  * @param qid 
//  * @param block 
//  * @param fid 
//  */
// void Ghost::PutGhost4Block_Wait(const qid_t* qid, GridBlock* block, const Field* fid) const {
//     //-------------------------------------------------------------------------
//     // reapply the boundary conditions from the latests updates I received
//     Compute4Block_Phys2Myself_(qid, block, fid);
//     //-------------------------------------------------------------------------
// }

// /**
//  * @brief copy the ghost values to my ghost area
//  */
// inline void Ghost::Compute4Block_Copy2Myself_(const ListGBLocal* ghost_list, const Field* fid, data_ptr data_trg) const {
//     //-------------------------------------------------------------------------
//     // the source block is always the same
//     const SubBlock block_src(M_GS, M_STRIDE, 0, M_N);

//     // loop on the ghost list
//     for (const auto gblock : (*ghost_list)) {
//         GridBlock* ngh_block = gblock->data_src();
//         // get info that change with the GP: where to put and from where to take it
//         MemLayout* block_trg = gblock;
//         data_ptr   data_src  = ngh_block->data(fid, ida_);
//         // copy the information
//         interp_->Copy(gblock->dlvl(), gblock->shift(), &block_src, data_src, block_trg, data_trg);
//     }
//     //-------------------------------------------------------------------------
// }

// /**
//  * @brief copy the ghost values to my coarse representation
//  */
// inline void Ghost::Compute4Block_Copy2Coarse_(const ListGBLocal* ghost_list, const Field* fid, mem_ptr ptr_trg) const {
//     //-------------------------------------------------------------------------
//     // the source block is always the same
//     const SubBlock block_src(M_GS, M_STRIDE, 0, M_N);

//     // the coarse block is computed later
//     lid_t coarse_start[3];
//     lid_t coarse_end[3];

//     // loop on the ghost list
//     for (const auto gblock : (*ghost_list)) {
//         GridBlock* ngh_block = gblock->data_src();
//         // set the coarse block to the correct position
//         for (lda_t id = 0; id < 3; id++) {
//             coarse_start[id] = interp_->CoarseFromBlock(gblock->start(id));
//             coarse_end[id]   = interp_->CoarseFromBlock(gblock->end(id));
//         }
//         SubBlock block_trg(interp_->CoarseNGhostFront(), interp_->CoarseStride(), coarse_start, coarse_end);
//         real_p   data_src = ngh_block->data(fid, ida_);
//         real_p   data_trg = ptr_trg + m_zeroidx(0, &block_trg);
//         // interpolate, the level is 1 coarser and the shift is unchanged
//         m_assert((gblock->dlvl() + 1) == 1 || (gblock->dlvl() + 1) == 0, "the difference of level MUST be 1");
//         interp_->Copy(gblock->dlvl() + 1, gblock->shift(), &block_src, data_src, &block_trg, data_trg);
//     }
//     //-------------------------------------------------------------------------
// }

// /**
//  * @brief get the ghost values to my ghost area using the @ref Wavelet::GetRma() function
//  */
// inline void Ghost::Compute4Block_GetRma2Myself_(const ListGBMirror* ghost_list, const Field* fid, data_ptr data_trg) const {
//     //-------------------------------------------------------------------------
//     // the source block is always the same
//     const SubBlock block_src(M_GS, M_STRIDE, 0, M_N);

//     // loop on the ghost list
//     for (const auto gblock : (*ghost_list)) {
//         MPI_Aint   disp_src  = gblock->data_src();
//         rank_t     disp_rank = gblock->rank();
//         MemLayout* block_trg = gblock;
//         // copy the information
//         interp_->GetRma(gblock->dlvl(), gblock->shift(), &block_src, disp_src, block_trg, data_trg, disp_rank, mirrors_window_);
//     }
//     //-------------------------------------------------------------------------
// }

// /**
//  * @brief get the ghost values to my coarse temp using the @ref Wavelet::GetRma() function
//  */
// inline void Ghost::Compute4Block_GetRma2Coarse_(const ListGBMirror* ghost_list, const Field* fid, mem_ptr ptr_trg) const {
//     m_assert(ghost_list != nullptr, "the ghost list cannot be null");
//     m_assert(ptr_trg != nullptr, "the ghost list cannot be null");
//     //-------------------------------------------------------------------------
//     // the source block is always the same
//     const SubBlock block_src(M_GS, M_STRIDE, 0, M_N);

//     // the coarse block is computed later
//     lid_t coarse_start[3];
//     lid_t coarse_end[3];

//     //................................................
//     // loop on the ghost list
//     for (const auto gblock : (*ghost_list)) {
//         MPI_Aint disp_src  = gblock->data_src();
//         rank_t   disp_rank = gblock->rank();
//         // set the coarse block to the correct position
//         for (lda_t id = 0; id < 3; id++) {
//             coarse_start[id] = interp_->CoarseFromBlock(gblock->start(id));
//             coarse_end[id]   = interp_->CoarseFromBlock(gblock->end(id));
//         }
//         SubBlock block_trg(interp_->CoarseNGhostFront(), interp_->CoarseStride(), coarse_start, coarse_end);
//         data_ptr data_trg = ptr_trg + m_zeroidx(0, &block_trg);

//         // interpolate, the level is 1 coarser and the shift is unchanged
//         m_assert((gblock->dlvl() + 1) == 1 || (gblock->dlvl() + 1) == 0, "the difference of level MUST be 1 or 0");
//         interp_->GetRma((gblock->dlvl() + 1), gblock->shift(), &block_src, disp_src, &block_trg, data_trg, disp_rank, mirrors_window_);
//     }
//     //-------------------------------------------------------------------------
// }

// /**
//  * @brief Apply the physical boundary conditions to my block
//  */
// inline void Ghost::Compute4Block_Phys2Myself_(const qid_t* qid, GridBlock* cur_block, const Field* fid) const {
//     //-------------------------------------------------------------------------
//     data_ptr data_trg = cur_block->data(fid, ida_);
//     for (auto gblock : (*cur_block->phys())) {
//         bctype_t bctype = fid->bctype(ida_, gblock->iface());
//         // get the correct face_start
//         if (bctype == M_BC_NEU) {
//             NeumanBoundary<M_WAVELET_N - 1> bc;
//             bc(gblock->iface(), face_start[gblock->iface()], cur_block->hgrid(), 0.0, gblock, data_trg);
//         } else if (bctype == M_BC_DIR) {
//             DirichletBoundary<M_WAVELET_N - 1> bc;
//             bc(gblock->iface(), face_start[gblock->iface()], cur_block->hgrid(), 0.0, gblock, data_trg);
//         } else if (bctype == M_BC_EXTRAP) {
//             ExtrapBoundary<M_WAVELET_N> bc;
//             bc(gblock->iface(), face_start[gblock->iface()], cur_block->hgrid(), 0.0, gblock, data_trg);
//         } else if (bctype == M_BC_ZERO) {
//             ZeroBoundary bc;
//             bc(gblock->iface(), face_start[gblock->iface()], cur_block->hgrid(), 0.0, gblock, data_trg);
//         } else {
//             m_assert(false, "this type of BC is not implemented yet or not valid %d", bctype);
//         }
//         //-------------------------------------------------------------------------
//     }
// }

// /**
//  * @brief interpolate my values to the coarse temp memory and apply the physical boundary conditions
//  */
// inline void Ghost::Compute4Block_Myself2Coarse_(const qid_t* qid, GridBlock* cur_block, const Field* fid, mem_ptr ptr_trg) const {
//     //-------------------------------------------------------------------------
//     SubBlock coarse_block(interp_->CoarseNGhostFront(), interp_->CoarseStride(), 0, M_HN);
//     // get memory details
//     const lid_t      shift[3]  = {0, 0, 0};
//     const MemLayout* block_src = cur_block;
//     const data_ptr   data_src  = cur_block->data(fid, ida_);
//     const MemLayout* block_trg = &coarse_block;
//     data_ptr         data_trg  = ptr_trg + m_zeroidx(0, &coarse_block);
//     // interpolate
//     interp_->Copy(1, shift, block_src, data_src, block_trg, data_trg);

//     //................................................
//     // do here some physics, to completely fill the coarse block before the interpolation
//     for (auto gblock : (*cur_block->phys())) {
//         // get the direction and the corresponding bctype
//         // const lda_t    dir    = gblock->dir();
//         const bctype_t bctype = fid->bctype(ida_, gblock->iface());
//         // in the face direction, the start and the end are already correct, only the fstart changes
//         lid_t fstart[3], coarse_start[3], coarse_end[3];
//         // in the other direction, we need to rescale the dimensions
//         for (lda_t id = 0; id < 3; id++) {
//             coarse_start[id] = interp_->CoarseFromBlock(gblock->start(id));
//             coarse_end[id]   = interp_->CoarseFromBlock(gblock->end(id));
//             fstart[id]       = interp_->CoarseFromBlock(face_start[gblock->iface()][id]);
//         }
//         // reset the coarse block and get the correct memory location
//         coarse_block.Reset(interp_->CoarseNGhostFront(), interp_->CoarseStride(), coarse_start, coarse_end);
//         data_ptr data_trg = ptr_trg + m_zeroidx(0, &coarse_block);
//         // get the correct face_start
//         if (bctype == M_BC_NEU) {
//             NeumanBoundary<M_WAVELET_N - 1> bc;
//             bc(gblock->iface(), fstart, cur_block->hgrid(), 0.0, &coarse_block, data_trg);
//         } else if (bctype == M_BC_DIR) {
//             DirichletBoundary<M_WAVELET_N - 1> bc;
//             bc(gblock->iface(), fstart, cur_block->hgrid(), 0.0, &coarse_block, data_trg);
//         } else if (bctype == M_BC_EXTRAP) {
//             ExtrapBoundary<M_WAVELET_N> bc;
//             bc(gblock->iface(), fstart, cur_block->hgrid(), 0.0, &coarse_block, data_trg);
//         } else if (bctype == M_BC_ZERO) {
//             ZeroBoundary bc;
//             bc(gblock->iface(), fstart, cur_block->hgrid(), 0.0, &coarse_block, data_trg);
//         } else {
//             m_assert(false, "this type of BC is not implemented yet or not valid %d", bctype);
//         }
//         //-------------------------------------------------------------------------
//     }
// }

// /**
//  * @brief refine from the coarse temp values to my ghost area, for local ghost blocks
//  */
// inline void Ghost::Compute4Block_Refine_(const ListGBLocal* ghost_list, const mem_ptr ptr_src, data_ptr data_trg) const {
//     //-------------------------------------------------------------------------
//     // the source block is always the same
//     const SubBlock block_src(interp_->CoarseNGhostFront(), interp_->CoarseStride(), -interp_->CoarseNGhostFront(), M_HN + interp_->CoarseNGhostBack());
//     const data_ptr data_src = ptr_src + m_zeroidx(0, &block_src);
//     lid_t          shift[3] = {0, 0, 0};

//     // loop on the ghost list
//     for (const auto gblock : (*ghost_list)) {
//         // interpolate, the level is 1 coarser and the shift is unchanged
//         interp_->Interpolate(-1, shift, &block_src, data_src, gblock, data_trg);
//     }
//     //-------------------------------------------------------------------------
// }

// /**
//  * @brief refine from the coarse temp values to my ghost area, for mirror ghost blocks
//  */
// inline void Ghost::Compute4Block_Refine_(const ListGBMirror* ghost_list, const mem_ptr ptr_src, data_ptr data_trg) const {
//     //-------------------------------------------------------------------------
//     // the source block is always the same
//     const SubBlock block_src(interp_->CoarseNGhostFront(), interp_->CoarseStride(), -interp_->CoarseNGhostFront(), M_HN + interp_->CoarseNGhostBack());
//     real_p         data_src = ptr_src + m_zeroidx(0, &block_src);
//     lid_t          shift[3] = {0, 0, 0};

//     // loop on the ghost list
//     for (const auto gblock : (*ghost_list)) {
//         // interpolate, the level is 1 coarser and the shift is unchanged
//         interp_->Interpolate(-1, shift, &block_src, data_src, gblock, data_trg);
//     }
//     //-------------------------------------------------------------------------
// }

// /**
//  * @brief coarsen to the coarse temp memory
//  */
// inline void Ghost::Compute4Block_Coarsen2Coarse_(const data_ptr data_src, mem_ptr ptr_trg) const {
//     //-------------------------------------------------------------------------
//     const lid_t shift[3] = {0, 0, 0};
//     // the target block is the full coarse tmp
//     const SubBlock block_trg(interp_->CoarseNGhostFront(), interp_->CoarseStride(), 0, M_HN);
//     data_ptr       data_trg = ptr_trg + m_zeroidx(0, &block_trg);
//     // the source block is the ghost extended block
//     const SubBlock block_src(M_GS, M_STRIDE, -interp_->nghost_front(), M_N + interp_->nghost_back());

//     // interpolate, the level is 1 coarser and the shift is unchanged
//     interp_->Interpolate(1, shift, &block_src, data_src, &block_trg, data_trg);
//     //-------------------------------------------------------------------------
// }

// /**
//  * @brief copy the data from the coarse to the parent's location
//  */
// inline void Ghost::Compute4Block_Copy2Parent_(const ListGBLocal* ghost_list, const mem_ptr ptr_src, const Field* fid) const {
//     //-------------------------------------------------------------------------
//     // the source block is always the same
//     const SubBlock block_src(interp_->CoarseNGhostFront(), interp_->CoarseStride(), 0, M_HN);
//     data_ptr       data_src = ptr_src + m_zeroidx(0, &block_src);

//     // loop on the ghost list
//     for (const auto gblock : (*ghost_list)) {
//         GridBlock* ngh_block = gblock->data_src();
//         data_ptr   data_trg  = ngh_block->data(fid, ida_);
//         // interpolate, the level is 1 coarser and the shift is unchanged
//         m_assert(gblock->dlvl() == 0, "we must have a level 0, here %d", gblock->dlvl());
//         interp_->Copy(gblock->dlvl(), gblock->shift(), &block_src, data_src, gblock, data_trg);
//     }
//     //-------------------------------------------------------------------------
// }

// /**
//  * @brief copy the data from the coarse to the parent location using RMA @ref Wavelet::PutRma()
//  */
// inline void Ghost::Compute4Block_PutRma2Parent_(const ListGBMirror* ghost_list, mem_ptr ptr_src) const {
//     //-------------------------------------------------------------------------
//     // the source block is always the same
//     const SubBlock block_src(interp_->CoarseNGhostFront(), interp_->CoarseStride(), 0, M_HN);

//     // loop on the ghost list
//     for (const auto gblock : (*ghost_list)) {
//         MPI_Aint disp_trg = gblock->data_src();
//         rank_t   trg_rank = gblock->rank();
//         // interpolate, the parent's mirror have been created to act on the tmp
//         m_assert(gblock->dlvl() == 0, "we must have a level 0, here %d", gblock->dlvl());
//         interp_->PutRma(gblock->dlvl(), gblock->shift(), &block_src, ptr_src, gblock, disp_trg, trg_rank, mirrors_window_);
//     }
//     //-------------------------------------------------------------------------
// }

/**
 * @brief Loop on the blocks that are mirrors and call a gop_t operation on them
 * 
 * @param op 
 * @param field 
 */
void Ghost::LoopOnMirrorBlock_(const gop_t op, const Field* field) {
    m_begin;
    m_assert(grid_->is_mesh_valid(), "mesh is not valid, unable to process");
    //-------------------------------------------------------------------------
    // get the grid info
    p8est_t*       forest = grid_->p4est_forest();
    p8est_ghost_t* ghost  = grid_->p4est_ghost();
    // const lid_t    nqlocal = ghost->mirrors.elem_count;  //number of ghost blocks

    //#pragma omp parallel for
    for (lid_t bid = 0; bid < n_mirror_to_send_; bid++) {
        // get the mirror quad, this is an empty quad (just a piggy3 struct)
        // p8est_quadrant_t* mirror = p8est_quadrant_array_index(&ghost->mirrors, local_to_mirrors[bid]);
        p8est_quadrant_t* mirror = p8est_quadrant_array_index(&ghost->mirrors, bid);
        p8est_tree_t*     tree   = p8est_tree_array_index(forest->trees, mirror->p.piggy3.which_tree);

        // build the mirror id
        qid_t myid;
        myid.cid = mirror->p.piggy3.local_num;   // cummulative id
        myid.mid = bid;                          // mirror id
        myid.tid = mirror->p.piggy3.which_tree;  // tree id
        // use it to retreive the actual quadrant in the correct tree
        p8est_quadrant_t* quad  = p8est_quadrant_array_index(&tree->quadrants, myid.cid - tree->quadrants_offset);
        GridBlock*        block = *(reinterpret_cast<GridBlock**>(quad->p.user_data));
        // send the task
        (this->*op)(&myid, block, field);
    }
    //-------------------------------------------------------------------------
    m_end;
}
