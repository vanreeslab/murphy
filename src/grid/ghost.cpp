#include "ghost.hpp"

#include <algorithm>
#include <limits>

#include "core/doop.hpp"
#include "core/macros.hpp"
#include "core/pointers.hpp"
#include "core/types.hpp"
#include "grid/boundary.hpp"
#include "omp.h"
#include "toolsp4est.hpp"
#include "wavelet/wavelet.hpp"

#define M_NNEIGHBOR 26

/**
 * @brief the localization of the interface
 */
// static lid_t face_start[6][3] = {{0, 0, 0}, {M_N, 0, 0}, {0, 0, 0}, {0, M_N, 0}, {0, 0, 0}, {0, 0, M_N}};

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
    m_verb("wavelet info:");
    m_verb("\t#ghost for refinement: %d %d", interp_->nghost_front_refine(), interp_->nghost_back_refine());
    m_verb("\t#ghost for coarsening: %d %d", interp_->nghost_front_coarsen(), interp_->nghost_back_coarsen());
    m_verb("\t#ghost for citerion and smoothing: %d %d", interp_->nghost_front_criterion_smooth(), interp->nghost_back_criterion_smooth());
    m_verb("\t#detail for criterion: %d %d", interp_->ndetail_citerion_extend_front(), interp_->ndetail_citerion_extend_back());
    m_verb("\t#detail for smoothing: %d %d", interp_->ndetail_smooth_extend_front(), interp_->ndetail_smooth_extend_back());
    m_verb("\tghost initialized with %s, nghost = %d %d, coarse nghost = %d %d", interp_->Identity().c_str(), interp_->nghost_front(), interp_->nghost_back(), interp_->CoarseNGhostFront(3), interp_->CoarseNGhostBack(3));

    // check that a fine block can provide enough ghosts to a coarse one
    m_assert(interp_->nghost_front() <= M_NHALF, "The memory for the ghost points is too small: M_NHALF = %d vs nghost = %d", M_NHALF, interp_->nghost_front());
    m_assert(interp_->nghost_back() <= M_NHALF, "The memory for the ghost points is too small: M_NHALF = %d vs nghost = %d", M_NHALF, interp_->nghost_back());
    // check that the N_GS is big enough
    m_assert(interp_->nghost_front() <= M_GS, "The memory for the ghost points is too small: M_NHALF = %d vs nghost = %d", M_NHALF, interp_->nghost_front());
    m_assert(interp_->nghost_back() <= M_GS, "The memory for the ghost points is too small: M_NHALF = %d vs nghost = %d", M_NHALF, interp_->nghost_back());

    //................................................
    // initialize the communications and the ghost's lists
    m_profStart(prof_, "ghost init comm");
    InitComm_();
    m_profStop(prof_, "ghost init comm");
    m_profStart(prof_, "ghost init list");
    InitList_();
    m_profStop(prof_, "ghost init list");

    //-------------------------------------------------------------------------
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
    m_profStart(prof_, "ghost free list");
    FreeList_();
    m_profStop(prof_, "ghost free list");

    m_profStart(prof_, "ghost free comm");
    FreeComm_();
    m_profStop(prof_, "ghost free comm");
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
    p8est_mesh_t*  mesh     = grid_->p4est_mesh();

    // sanity checks
    m_assert(mpi_comm == MPI_COMM_WORLD, "the comm should be a comm world");

    //................................................
    // allocate the local2disp and the status window
    MPI_Info info;
    MPI_Info_create(&info);
    MPI_Info_set(info, "no_locks", "true");

    // displacement
    m_profStart(prof_, "MPI_Win_create");
    MPI_Aint  win_disp_mem_size = mesh->local_num_quadrants * sizeof(MPI_Aint);
    MPI_Aint* local2disp        = static_cast<MPI_Aint*>(m_calloc(win_disp_mem_size));
    MPI_Win   local2disp_window = MPI_WIN_NULL;
    MPI_Win_create(local2disp, win_disp_mem_size, sizeof(MPI_Aint), info, MPI_COMM_WORLD, &local2disp_window);
    m_assert(win_disp_mem_size >= 0, "the memory size should be >=0");
    m_assert(local2disp_window != MPI_WIN_NULL, "window must be ready");
    m_verb("allocating %ld bytes in the window for %d active quad", win_disp_mem_size, mesh->local_num_quadrants);
    m_profStop(prof_, "MPI_Win_create");

    // free the info
    MPI_Info_free(&info);

    //................................................
    // compute the number of admissible local mirrors and store their reference in the array
    iblock_t    active_mirror_count = 0;
    const lid_t nmlocal             = ghost->mirrors.elem_count;  //number of ghost blocks
    // get the base address
    // fill the displacement
    for (iblock_t im = 0; im < nmlocal; im++) {
        qdrt_t* mirror       = p8est_quadrant_array_index(&ghost->mirrors, im);
        level_t mirror_level = p4est_GetQuadFromMirror(forest, mirror)->level;
        m_assert(mirror->level == mirror_level, "the two levels must be the same: %d vs %d", mirror->level, mirror_level);
        // update the counters if the mirror is admissible
        // we need to have called the function p4est_balance()!!
        if ((min_level_ - 1) <= mirror_level && mirror_level <= (max_level_ + 1)) {
            // store the displacement in the local2disp_ array
            iblock_t local_id = mirror->p.piggy3.local_num;
            m_assert(local_id >= 0, "the address cannot be negative");
            m_assert(local_id < mesh->local_num_quadrants, "the address cannot be > %d", mesh->local_num_quadrants);
            local2disp[local_id] = active_mirror_count * CartBlockMemNum(1);
            active_mirror_count++;
            m_assert(active_mirror_count <= nmlocal, "the number of mirrors cannot be bigger than the local number:  %d vs %d", active_mirror_count, nmlocal);
        }
    }

    //................................................
    // post the exposure epoch and start the access one for local2mirrors
    m_assert(ingroup_ != MPI_GROUP_NULL, "call the InitComm function first!");
    m_assert(outgroup_ != MPI_GROUP_NULL, "call the InitComm function first!");

    // start the exposure epochs if any
    m_profStart(prof_, "init list on blocks");
    MPI_Win_post(ingroup_, 0, local2disp_window);
    MPI_Win_start(outgroup_, 0, local2disp_window);

    // init the list on every active block that matches the level requirements
    const ForestGrid* mygrid = grid_;
    for (level_t il = min_level_; il <= max_level_; il++) {
        DoOpMeshLevel(nullptr, &GridBlock::GhostInitLists, grid_, il, mygrid, interp_, local2disp_window);
    }

    // complete the epoch and wait for the exposure one
    MPI_Win_complete(local2disp_window);
    MPI_Win_wait(local2disp_window);
    m_profStop(prof_, "init list on blocks");

    //................................................
    m_profStart(prof_,"MPI_Win_free");
    MPI_Win_free(&local2disp_window);
    m_free(local2disp);
    m_profStop(prof_,"MPI_Win_free");

    //................................................
    // allocate the status

    //-------------------------------------------------------------------------
    m_verb("Ghost lists initialization is done");
    m_end;
}

void Ghost::FreeList_() {
    m_begin;
    //-------------------------------------------------------------------------
    for (level_t il = min_level_; il <= max_level_; il++) {
        DoOpMeshLevel(nullptr, &GridBlock::GhostFreeLists, grid_, il);
    }
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief Initialize the communication pattern: init the windows for RMA access
 * 
 * We initiate the window mirrors_window_ and the associated groups: ingroup_ & outgroup_
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
    p8est_mesh_t*  mesh     = grid_->p4est_mesh();

    //................................................
    // compute the number of admissible local mirrors and store their reference in the array
    iblock_t n_mirror_to_send = 0;
    for (iblock_t im = 0; im < ghost->mirrors.elem_count; im++) {
        qdrt_t* mirror       = p8est_quadrant_array_index(&ghost->mirrors, im);
        level_t mirror_level = p4est_GetQuadFromMirror(forest, mirror)->level;
        // update the counters if the mirror is admissible (i.e. it satisfies the requirements)
        // we need to have called the function p4est_balance()!!
        if ((min_level_ - 1) <= mirror_level && mirror_level <= (max_level_ + 1)) {
            n_mirror_to_send++;
        }
    }
    m_verb("I have %d mirrors to send", n_mirror_to_send);

    // allocate memory
    m_profStart(prof_, "allocate mem");
    MPI_Aint win_mem_size = n_mirror_to_send * CartBlockMemNum(1) * sizeof(real_t);
    mirrors_              = static_cast<real_t*>(m_calloc(win_mem_size));
    MPI_Aint win_status_mem_size = mesh->local_num_quadrants * sizeof(short_t);
    status_                      = static_cast<short_t*>(m_calloc(win_status_mem_size));
    m_profStop(prof_, "allocate mem");

    // initialize the Window by allocating the memory space needed for the exchange
    m_profStart(prof_, "MPI_Win_create");
    MPI_Info info;
    MPI_Info_create(&info);
    MPI_Info_set(info, "no_locks", "true");

    // window for the mirrors
    // because of alignement issues (https://github.com/open-mpi/ompi/issues/7955), cannot use this one
    // MPI_Win_allocate(win_mem_size, sizeof(real_t), info, mpi_comm, &mirrors_, &mirrors_window_);
    MPI_Win_create(mirrors_, win_mem_size, sizeof(real_t), info, MPI_COMM_WORLD, &mirrors_window_);
    MPI_Win_create(status_, win_status_mem_size, sizeof(short_t), info, MPI_COMM_WORLD, &status_window_);

    MPI_Info_free(&info);
    m_profStop(prof_, "MPI_Win_create");

    // get an idea of the percentage of blocks
#ifndef NDEBUG
    {
        real_t   global_ratio = 0.0;
        real_t   ratio        = 0.0;
        iblock_t n_quad_local = mesh->local_num_quadrants;
        if (n_quad_local > 0) {
            ratio = static_cast<real_t>(n_mirror_to_send) / static_cast<real_t>(n_quad_local);
        }
        MPI_Allreduce(&ratio, &global_ratio, 1, M_MPI_REAL, MPI_SUM, MPI_COMM_WORLD);
        int comm_size = 0;
        MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
        m_log("mirrors are %2.2f %% of the blocks", 100.0 * global_ratio / (comm_size));

        m_assert(m_isaligned(mirrors_), "the mirror temp array is not aligned");
        m_assert(mirrors_window_ != MPI_WIN_NULL, "the MPI window created is null, which is not a good news");
        m_assert(status_window_ != MPI_WIN_NULL, "the MPI window created is null, which is not a good news");
    }
#endif

    //................................................
    // get the list of ranks that will generate a call to access my mirrors
    rank_t n_in_group = 0;
    // over allocate the array to its max size and fill it only partially... not great
    rank_t* group_ranks = static_cast<rank_t*>(m_calloc(mpi_size * sizeof(rank_t)));
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
                group_ranks[n_in_group] = ir;
                n_in_group += 1;
                break;
            }
        }
    }
    // get the RMA mirror group ready - the group that will need my info
    MPI_Group global_group;
    MPI_Comm_group(MPI_COMM_WORLD, &global_group);
    // MPI_Win_get_group(mirrors_window_, &win_group);
    MPI_Group_incl(global_group, n_in_group, group_ranks, &ingroup_);
    m_assert(!(n_in_group == 0 && ingroup_ != MPI_GROUP_EMPTY), "if there is no cpu, the group must be empty");

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
                group_ranks[n_in_group] = ir;
                n_in_group += 1;
                break;
            }
        }
    }
    // add the cpus that will get a call from me
    MPI_Group_incl(global_group, n_in_group, group_ranks, &outgroup_);
    MPI_Group_free(&global_group);
    m_assert(!(n_in_group == 0 && ingroup_ != MPI_GROUP_EMPTY), "if there is no cpu, the group must be empty");
#ifndef NDEBUG
    {
        int test;
        MPI_Group_compare(ingroup_, outgroup_, &test);
        m_assert(test != MPI_UNEQUAL, "the ingroup and outgroup must be the same: test = %d vs %d, %d, %d and %d", test, MPI_IDENT, MPI_CONGRUENT, MPI_SIMILAR, MPI_UNEQUAL);
    }
#endif

    //................................................
    // free the allocated memory
    m_free(group_ranks);
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief Free the communication window mirrors_window_ and the groups ingroup_ & outgroup_
 * 
 */
void Ghost::FreeComm_() {
    m_begin;
    //-------------------------------------------------------------------------
    // free the group
    m_profStart(prof_, "MPI_Group_free");
    MPI_Group_free(&ingroup_);
    MPI_Group_free(&outgroup_);
    m_profStop(prof_, "MPI_Group_free");
    // free the window
    m_profStart(prof_, "MPI_Win_free");
    MPI_Win_free(&mirrors_window_);
    MPI_Win_free(&status_window_);
    m_profStop(prof_, "MPI_Win_free");

    m_profStart(prof_, "free mem");
    m_free(mirrors_);
    m_free(status_);
    m_profStop(prof_, "free mem");
    //-------------------------------------------------------------------------
    m_end;
}

void Ghost::UpdateStatus() {
    m_begin;
    //-------------------------------------------------------------------------
    m_profStart(prof_, "update status");
    m_profStart(prof_, "fill status");
    // get the status to the array
    for (level_t il = min_level_; il <= max_level_; il++) {
        DoOpMeshLevel(nullptr, &GridBlock::SetNewByCoarsening, grid_, il, status_);
    }
    m_profStop(prof_, "fill status");

    // start the exposure epochs if any (we need to be accessed by the neighbors even is we have not block on that level)
    m_profStart(prof_, "comm PS");
    // if (ingroup_ != MPI_GROUP_EMPTY) {
    MPI_Win_post(ingroup_, 0, status_window_);
    // }
    // if (outgroup_ != MPI_GROUP_EMPTY) {
    MPI_Win_start(outgroup_, 0, status_window_);
    // }
    m_profStop(prof_, "comm PS");

    m_profStart(prof_, "RMA Get");
    // update neigbbor status, only use the already computed status on level il + 1
    for (level_t il = min_level_; il <= max_level_; il++) {
        DoOpMeshLevel(nullptr, &GridBlock::GetNewByCoarseningFromNeighbors, grid_, il, status_, status_window_);
    }
    m_profStop(prof_, "RMA Get");

    // close the access epochs
    m_profStart(prof_, "comm CW");
    // if (outgroup_ != MPI_GROUP_EMPTY) {
    MPI_Win_complete(status_window_);
    // }
    // if (ingroup_ != MPI_GROUP_EMPTY) {
    MPI_Win_wait(status_window_);
    // }
    m_profStop(prof_, "comm CW");
    m_profStop(prof_, "update status");
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief sets the length of ghosting for all the blocks and adapt the required length if needed
 * 
 * We must be able to overwrite and coarsen while doing the ghosting if the grid contains different levels.
 * Therefore we have to guarantee that the have enough ghosts points.
 * 
 * @param ghost_len the desired length of ghosts, might be changed after the function returns
 */
void Ghost::SetLength(bidx_t ghost_len[2]) {
    m_profStart(prof_, "set length");
    //-------------------------------------------------------------------------
    // adapt the ghost lengths if we are a MR grid

    m_profStart(prof_, "loop");
    const bool   is_grid_mr = grid_->MaxLevel() > grid_->MinLevel();
    const bidx_t new_len[2] = {m_max(ghost_len[0], is_grid_mr * m_max(interp_->nghost_front_overwrite(), interp_->nghost_front_coarsen())),
                               m_max(ghost_len[1], is_grid_mr * m_max(interp_->nghost_back_overwrite(), interp_->nghost_back_coarsen()))};
    m_profStop(prof_, "loop");
    m_assert(new_len[0] <= M_GS, "there is not enough space for the ghosts -> requested: %d, actual: %d, space; %d", ghost_len[0], new_len[0], M_GS);
    m_assert(new_len[1] <= M_GS, "there is not enough space for the ghosts -> requested: %d, actual: %d, space; %d", ghost_len[1], new_len[1], M_GS);
    m_verb("length set: %s %d %d", (new_len[0] == ghost_len[0] && new_len[1] == ghost_len[1]) ? "" : "!WARNING! the ghost lenghts have been changed to", new_len[0], new_len[1]);
    
    m_profStart(prof_, "loop");
    for (level_t il = min_level_; il <= max_level_; il++) {
        DoOpMeshLevel(nullptr, &GridBlock::GhostUpdateSize, grid_, il, new_len);
    }
    m_profStop(prof_, "loop");

    // overwrites the real lengths
    ghost_len[0] = new_len[0];
    ghost_len[1] = new_len[1];

    // just store a small string for the prof
    prof_msg_ = " (" + std::to_string(ghost_len[0]) + "," + std::to_string(ghost_len[1]) + ")";
    //-------------------------------------------------------------------------
    m_profStop(prof_, "set length");
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
    cur_ida_ = ida;
    m_profStart(prof_, "pullghost post" + prof_msg_);
    //................................................
    // fill the Window memory with the Mirror information
    m_profStart(prof_, "(01) push to window");
    LoopOnMirrorBlock_(&Ghost::PushToWindow4Block, field);
    m_profStop(prof_, "(01) push to window");

    //................................................
    // post the exposure epoch for my own mirrors: I am a target warning that origin group will RMA me
    // start the access epoch, to get info from neighbors: I am an origin warning that I will RMA the target group
    m_profStart(prof_, "(02) RMA - post start");
    MPI_Win_post(ingroup_, 0, mirrors_window_);
    MPI_Win_start(outgroup_, 0, mirrors_window_);
    m_profStop(prof_, "(02) RMA - post start");

    m_profStart(prof_, "(03) RMA - get");
    for (level_t il = min_level_; il <= max_level_; il++) {
        DoOpMeshLevel(nullptr, &GridBlock::GhostGet_Post, grid_, il, field, ida, interp_, mirrors_window_);
    }
    m_profStop(prof_, "(03) RMA - get");

    //................................................
    // start what can be done = sibling and parents local copy + physical BC + myself copy
    m_profStart(prof_, "(04) computation");
    for (level_t il = min_level_; il <= max_level_; il++) {
        DoOpMeshLevel(nullptr, &GridBlock::GhostGet_Cmpt, grid_, il, field, ida, interp_);
    }
    m_profStop(prof_, "(04) computation");
    m_profStop(prof_, "pullghost post" + prof_msg_);
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
    m_assert(cur_ida_ == ida, "the ongoing dimension (%d) must be over first", cur_ida_);
    m_assert(grid_->is_mesh_valid(), "the mesh needs to be valid before entering here");
    //-------------------------------------------------------------------------
    m_profStart(prof_, "pullghost wait" + prof_msg_);
    //................................................
    // finish the access epochs for the exposure epoch to be over
    m_profStart(prof_, "(05) RMA - complete");
    MPI_Win_complete(mirrors_window_);
    m_profStop(prof_, "(05) RMA - complete");
    m_profStart(prof_, "(06) RMA - wait");
    MPI_Win_wait(mirrors_window_);
    m_profStop(prof_, "(06) RMA - wait");

    // we now have all the information needed to compute the ghost points in coarser blocks
    m_profStart(prof_, "(07) computation");
    for (level_t il = min_level_; il <= max_level_; il++) {
        DoOpMeshLevel(nullptr, &GridBlock::GhostGet_Wait, grid_, il, field, ida, interp_);
    }
    m_profStop(prof_, "(07) computation");

    //................................................

    // post exposure and access epochs for to put the values to my neighbors
    m_profStart(prof_, "(08) RMA - post start");
    MPI_Win_post(ingroup_, 0, mirrors_window_);
    MPI_Win_start(outgroup_, 0, mirrors_window_);
    m_profStop(prof_, "(08) RMA - post start");

    // start what can be done = sibling and parents copy
    m_profStart(prof_, "(09) RMA - put");
    for (level_t il = min_level_; il <= max_level_; il++) {
        DoOpMeshLevel(nullptr, &GridBlock::GhostPut_Post, grid_, il, field, ida, interp_, mirrors_window_);
    }
    m_profStop(prof_, "(09) RMA - put");

    // finish the access epochs for the exposure epoch to be over
    m_profStart(prof_, "(10) RMA complete");
    MPI_Win_complete(mirrors_window_);
    m_profStop(prof_, "(10) RMA complete");
    m_profStart(prof_, "(11) RMA wait");
    MPI_Win_wait(mirrors_window_);
    m_profStop(prof_, "(11) RMA wait");

    m_profStart(prof_, "(12) pull from window");
    LoopOnMirrorBlock_(&Ghost::PullFromWindow4Block, field);
    m_profStop(prof_, "(12) pull from window");

    // we now have all the information needed, we finish with a physbc
    m_profStart(prof_, "(13) computation");
    for (level_t il = min_level_; il <= max_level_; il++) {
        DoOpMeshLevel(nullptr, &GridBlock::GhostPut_Wait, grid_, il, field, ida, interp_);
    }
    m_profStop(prof_, "(13) computation");
    m_profStop(prof_, "pullghost wait" + prof_msg_);
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief Push the local mirrors to the window's associated memory
 * 
 * @param qid the quarant id considered
 * @param block the grid block considered
 * @param fid the field ID
 */
void Ghost::PushToWindow4Block(const qid_t* qid, GridBlock* block, const Field* fid) const {
    m_assert(cur_ida_ >= 0, "the current working dimension has to be correct");
    m_assert(cur_ida_ < fid->lda(), "the current working dimension has to be correct");
    //-------------------------------------------------------------------------
    // recover the mirro spot using the mirror id
    real_p  mirror = mirrors_ + qid->mid * CartBlockMemNum(1);
    mem_ptr data   = block->pointer(fid, cur_ida_);
    // m_assume_aligned(mirror);
    // m_assume_aligned(data);
    memcpy(mirror, data(), CartBlockMemNum(1) * sizeof(real_t));
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
    m_assert(cur_ida_ >= 0, "the current working dimension has to be correct");
    m_assert(cur_ida_ < fid->lda(), "the current working dimension has to be correct");
    //-------------------------------------------------------------------------
    real_p   mirror = mirrors_ + qid->mid * CartBlockMemNum(1) + m_zeroidx(0, block);
    data_ptr data   = block->data(fid, cur_ida_);
    // m_assume_aligned(mirror);
    // m_assume_aligned(data);

    for (auto* gblock : (*block->ghost_children())) {
        const lid_t start[3] = {gblock->start(0), gblock->start(1), gblock->start(2)};
        const lid_t end[3]   = {gblock->end(0), gblock->end(1), gblock->end(2)};

        real_t* data_src = mirror + m_idx(start[0], start[1], start[2], 0, block->stride());
        real_t* data_trg = data.Write(start[0], start[1], start[2], 0, block);

        // copy the value = sendrecv to myself to the correct spot
        MPI_Status   status;
        MPI_Datatype dtype;
        ToMPIDatatype(start, end, block->stride(), 1, &dtype);
        MPI_Sendrecv(data_src, 1, dtype, 0, 0, data_trg, 1, dtype, 0, 0, MPI_COMM_SELF, &status);
        MPI_Type_free(&dtype);
    }
    //-------------------------------------------------------------------------
}

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
    p8est_t*       forest  = grid_->p4est_forest();
    p8est_ghost_t* ghost   = grid_->p4est_ghost();
    const lid_t    nqlocal = ghost->mirrors.elem_count;  //number of ghost blocks

    //#pragma omp parallel for
    for (lid_t bid = 0; bid < nqlocal; bid++) {
        // get the mirror quad, this is an empty quad (just a piggy3 struct)
        // p8est_quadrant_t* mirror = p8est_quadrant_array_index(&ghost->mirrors, local_to_mirrors[bid]);
        p8est_quadrant_t* mirror = p8est_quadrant_array_index(&ghost->mirrors, bid);
        p8est_tree_t*     tree   = p8est_tree_array_index(forest->trees, mirror->p.piggy3.which_tree);

        if ((min_level_ - 1) <= mirror->level && mirror->level <= (max_level_ + 1)) {
            // build the mirror id
            qid_t myid;
            myid.cid = mirror->p.piggy3.local_num;   // cummulative id
            myid.mid = bid;                          // mirror id
            myid.tid = mirror->p.piggy3.which_tree;  // tree id
            // use it to retreive the actual quadrant in the correct tree
            p8est_quadrant_t* quad = p8est_quadrant_array_index(&tree->quadrants, myid.cid - tree->quadrants_offset);
            // GridBlock*        block = *(reinterpret_cast<GridBlock**>(quad->p.user_data));
            GridBlock* block = p4est_GetGridBlock(quad);
            // send the task
            (this->*op)(&myid, block, field);
        }
    }
    //-------------------------------------------------------------------------
    m_end;
}
