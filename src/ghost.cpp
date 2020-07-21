#include "ghost.hpp"

#include <limits>

#include "boundary.hpp"
#include "murphy.hpp"
#include "omp.h"
#include "toolsp4est.hpp"
#include "wavelet.hpp"


#define M_NNEIGHBOR 26

/**
 * @brief returns the number of ghost points for the coarse block (in the front)
 */
static lid_t CoarseNGhostFront(Interpolator* interp) {
    // we need the max between the number of scaling in front 
    //      = (interp->nghost_front() / 2)
    // and the number of details
    //      = (interp->nghost_front() / 2) + interp->nrefine_front()
    const lid_t gp = (interp->nghost_front() / 2) + interp->nrefine_front();
    return gp;
}

/**
 * @brief returns the number of ghost points for the coarse block (in the back)
 */
static lid_t CoarseNGhostBack(Interpolator* interp) {
    // we need the max between the number of scaling
    //      = (interp->nghost_back()+1) / 2
    const lid_t gp_scaling = ((interp->nghost_back()+1) / 2);
    // and the number of details
    //      = (interp->nghost_back()) / 2 +  + interp->nrefine_back()
    const lid_t gp_detail  = ((interp->nghost_back()) / 2) + interp->nrefine_back();
    return m_max(gp_scaling,gp_detail);
}

/**
 * @brief returns the stride of the coarse block
 * we need ghost points in front, at the back and M_HN points inbetween
 */
static size_t CoarseStride(Interpolator* interp) {
    const lid_t gp_front = CoarseNGhostFront(interp);
    const lid_t gp_back  = CoarseNGhostBack(interp);
    // the stride = ghost point in front + in back + M_HN in between
    return gp_front + M_HN + gp_back;
}

/**
 * @brief return the memory size (in byte!) of the coarse block
 * 
 * @param interp 
 * @return size_t 
 */
static size_t CoarseMemSize(Interpolator* interp) {
    return CoarseStride(interp) * CoarseStride(interp) * CoarseStride(interp) * sizeof(real_t);
}

/**
 * @brief Given a face, edge or a corner returns the outgoing normal (sign)
 * 
 * @param ibidule for a face (`ibidule<6`), an edge (`6<= ibidule < 18`) or a corner (`18<= ibidule < 26`)
 * @param sign the sign of the outgoing normal
 */
static void GhostGetSign(sid_t ibidule, real_t sign[3]) {
    // we need to find the sign = the direction of the normal:
    sign[0] = 0.0;
    sign[1] = 0.0;
    sign[2] = 0.0;

    // check depending on the plane, the edge of the corner
    if (ibidule < 6) {
        sid_t dir = ibidule / 2;
        sign[dir] = ((ibidule % 2) == 1) ? 1.0 : -1.0;
    } else if (ibidule < 18) {
        sid_t iedge = ibidule - 6;
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
        sid_t dir  = iedge / 4;           // this is the direction of the edge
        sid_t dir1 = (dir == 0) ? 1 : 0;  // dir1 in the plane: dir1 = x if dir = y or z or y if dir = x
        sid_t dir2 = (dir == 2) ? 1 : 2;  // dir2 in the plane: dir2 = y if dir=z, = z if dir=x or dir = y
        sign[dir1] = ((iedge % 4) % 2) == 1 ? +1.0 : -1.0;
        sign[dir2] = ((iedge % 4) / 2) == 1 ? +1.0 : -1.0;
    } else {
        sid_t icorner = ibidule - 18;
        sign[0]       = (icorner % 2) == 1 ? +1.0 : -1.0;
        sign[1]       = ((icorner % 4) / 2) == 1 ? +1.0 : -1.0;
        sign[2]       = (icorner / 4) == 1 ? +1.0 : -1.0;
    }

    m_assert(sign[0] == 0.0 || sign[0] == 1.0 || sign[0] == -1.0, "wrong sign value: %e", sign[0]);
    m_assert(sign[1] == 0.0 || sign[1] == 1.0 || sign[1] == -1.0, "wrong sign value: %e", sign[1]);
    m_assert(sign[2] == 0.0 || sign[2] == 1.0 || sign[2] == -1.0, "wrong sign value: %e", sign[2]);
}

/**
 * @brief given a starting a ghost range for a block, transform it into the needed range for its coarse representation
 * if the point is inside the block, we divide its id by 2, if the point is in the ghost points,
 * we take ALL the coarse ghost points
 * 
 * 
 * @warning this is magic...
 * 
 * we compute:
 * b = a + M_N, such that
 *      if a is in the ghost points, b < M_N
 *      if a is in the center points, M_N <= b < 2* M_N
 *      if a is in the ghost points, 2* M_N <= b
 * 
 * c is 0,1,2,3:
 *      = 0 if a is in the negative ghost points
 *      = 1 if a is in the center points, including 0
 *      = 2 if a is M_N
 *      = 3 if a is in the negative GP
 * 
 * the correct ID is returned based on the value of c
 * 
 * @param a the id in the block
 * @param interp the interpolator used
 * @return lid_t 
 */
static inline lid_t CoarseFromBlock(const lid_t a, Interpolator* interp) {
    const lid_t gp_front = CoarseNGhostFront(interp);
    const lid_t gp_back  = CoarseNGhostBack(interp);
    const lid_t b        = (a + M_N);
    const lid_t c        = (b / M_N) + (a > M_N);
    const lid_t res[4]   = {-gp_front, (a / 2), M_HN, M_HN + gp_back};
    // return the correct choice
    return res[c];
}

void CallGhostInitList(const qid_t* qid, GridBlock* block, nullptr_t fid, Ghost* ghost) {
    m_assert(fid == nullptr, "this pointer has to be null");
    ghost->InitList4Block(qid, block);
}
void CallGetGhost4Block_Post(const qid_t* qid, GridBlock* block, Field* fid, Ghost* ghost) {
    ghost->GetGhost4Block_Post(qid, block, fid);
}
void CallGetGhost4Block_Wait(const qid_t* qid, GridBlock* block, Field* fid, Ghost* ghost) {
    ghost->GetGhost4Block_Wait(qid, block, fid);
}
void CallPutGhost4Block_Post(const qid_t* qid, GridBlock* block, Field* fid, Ghost* ghost) {
    ghost->PutGhost4Block_Post(qid, block, fid);
}
void CallPutGhost4Block_Wait(const qid_t* qid, GridBlock* block, Field* fid, Ghost* ghost) {
    ghost->PutGhost4Block_Wait(qid, block, fid);
}
// void CallGhostPullFromGhost(const qid_t* qid, GridBlock* block, Field* fid, Ghost* ghost) {
//     //-------------------------------------------------------------------------
//     ghost->PullFromGhost4Block(qid, block, fid);
//     //-------------------------------------------------------------------------
// }



/**
 * @brief Construct a new Ghost object 
 * 
 * see Ghost::Ghost(ForestGrid* grid, const level_t min_level, const level_t max_level, Interpolator* interp) for details
 * 
 * @param grid the ForestGrid to use, must have been initiated using ForestGrid::SetupP4estGhostMesh() 
 * @param interp the interpolator to use, will drive the number of ghost points to consider
 */
Ghost::Ghost(ForestGrid* grid, Interpolator* interp) : Ghost(grid, -1, P8EST_MAXLEVEL + 1, interp) {
    // simply call the detailed constructor
}

/**
 * @brief Construct a new Ghost, allocate the ghost lists and initiates the communications
 * 
 * Once created, the ghost is fixed for a given grid. if the grid changes, a new Ghost objects has to be created.
 * 
 * @note: The Interpolator has to be given beforehands because the it drives the number of actual GP to consider.
 * While for memory alignement, the number of ghost points is given by M_GS, the wavelet does not require that many ghost points to
 * be computed. To reduce the memory cost, only the needed ghost points will be computed.
 * 
 * @param grid the ForestGrid to use, must have been initiated using ForestGrid::SetupP4estGhostMesh() 
 * @param min_level the minimum level on which the GP are initiated
 * @param max_level the maximum level on which the GP are initiated
 * @param interp the interpolator to use, will drive the number of ghost points to consider
 */
Ghost::Ghost(ForestGrid* grid, const level_t min_level, const level_t max_level, Interpolator* interp) {
    m_begin;
    m_assert(grid->is_mesh_valid(), "the mesh needs to be valid before entering here");
    //-------------------------------------------------------------------------
    // get the important pointers
    grid_              = grid;
    interp_            = interp;
    // store the level information
    min_level_ = m_max(min_level, 0);
    max_level_ = m_min(max_level, P8EST_QMAXLEVEL);

    //................................................
    // get how many active quads should be considered and allocate the ghost ptr
    p8est_mesh_t* mesh = grid_->mesh();
    n_active_quad_     = 0;
    for (level_t il = min_level_; il <= max_level_; il++) {
        n_active_quad_ += p4est_NumQuadOnLevel(mesh, il);
    }

    // store the number of ghosts needed
    nghost_[0] = interp_->nghost_front();
    nghost_[1] = interp_->nghost_back();
    m_assert(nghost_[0] <= M_GS, "The memory for the ghost points is too small: M_GS = %d vs nghost = %d", M_GS, nghost_[0]);
    m_assert(nghost_[1] <= M_GS, "The memory for the ghost points is too small: M_GS = %d vs nghost = %d", M_GS, nghost_[1]);

    //................................................
    // initialize the communications and the ghost's lists
    InitComm_();
    InitList_();

    //-------------------------------------------------------------------------
    m_log("ghost initialized with %s, nghost = %d %d, coarse ghost = %d %d", interp_->Identity().c_str(), nghost_[0], nghost_[1], CoarseNGhostFront(interp_), CoarseNGhostBack(interp_));
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

    // // end the requests
    // for (int ir = 0; ir < n_send_request_; ir++) {
    //     MPI_Request_free(mirror_send_ + ir);
    // }
    // m_free(mirror_send_);
    // for (int ir = 0; ir < n_recv_request_; ir++) {
    //     MPI_Request_free(ghost_recv_ + ir);
    // }
    // m_free(ghost_recv_);

    // // free the memory bylevel
    // m_free(local_to_mirrors);
    // m_free(ghost_to_local_);

    // // free the temp memory
    // int nthreads = omp_get_max_threads();
    // for (int it = 0; it < nthreads; it++) {
    //     m_free(coarse_tmp_[it]);
    // }
    // m_free(coarse_tmp_);

    // // free the data memory
    // m_free(mirrors_);
    // m_free(ghosts_);
    //-------------------------------------------------------------------------
    m_end;
}


void Ghost::InitList_(){
    m_begin;
    //-------------------------------------------------------------------------
    // allocate the lists for the corresponding quads
    block_children_ = (ListGBLocal**)m_calloc(n_active_quad_ * sizeof(ListGBLocal*));
    block_sibling_  = (ListGBLocal**)m_calloc(n_active_quad_ * sizeof(ListGBLocal*));
    block_parent_   = (ListGBLocal**)m_calloc(n_active_quad_ * sizeof(ListGBLocal*));
    ghost_children_ = (ListGBMirror**)m_calloc(n_active_quad_ * sizeof(ListGBMirror*));
    ghost_sibling_  = (ListGBMirror**)m_calloc(n_active_quad_ * sizeof(ListGBMirror*));
    ghost_parent_   = (ListGBMirror**)m_calloc(n_active_quad_ * sizeof(ListGBMirror*));
    phys_           = (listGBPhysic**)m_calloc(n_active_quad_ * sizeof(listGBPhysic*));
    // init the lists
    for (int ib = 0; ib < n_active_quad_; ib++) {
        // purge everything
        block_children_[ib] = new ListGBLocal();
        block_sibling_[ib]  = new ListGBLocal();
        block_parent_[ib]   = new ListGBLocal();
        ghost_children_[ib] = new ListGBMirror();
        ghost_sibling_[ib]  = new ListGBMirror();
        ghost_parent_[ib]   = new ListGBMirror();
        phys_[ib]           = new listGBPhysic();
    }

    //................................................
    // get stupid MPI info
    int            mpi_size = grid_->mpisize();
    MPI_Comm       mpi_comm = grid_->mpicomm();
    p8est_t*       forest   = grid_->forest();
    p8est_ghost_t* ghost    = grid_->ghost();

    // allocate the array to link the local_id to the mirror displacement and init it
    MPI_Info info;
    MPI_Info_create(&info);
    MPI_Info_set(info, "no_locks", "true");
    MPI_Aint win_mem_size = n_active_quad_ * sizeof(MPI_Aint);
    m_log("allocating %ld bytes in the window for %d active quad",win_mem_size,n_active_quad_);
    MPI_Win_allocate(win_mem_size, sizeof(MPI_Aint), info, mpi_comm, &local2disp_, &local2disp_window_);
    MPI_Info_free(&info);
    for (iblock_t ib = 0; ib < n_active_quad_; ib++) {
        local2disp_[ib] = 0;
    }
    MPI_Win_fence(0, local2disp_window_);

    //................................................
    // compute the number of admissible local mirros and store their reference in the array
    iblock_t count = 0;
    for (iblock_t im = 0; im < ghost->mirrors.elem_count; im++) {
        qdrt_t* mirror       = p8est_quadrant_array_index(&ghost->mirrors, im);
        sid_t   mirror_level = p4est_GetQuadFromMirror(forest, mirror)->level;
        // update the counters if the mirror is admissible
        // we need to have called the function p4est_balance()!!
        if ((min_level_ - 1) <= mirror_level && mirror_level <= (max_level_ + 1)) {
            // store the displacement in the local2disp_ array
            iblock_t local_id     = mirror->p.piggy3.local_num;
            local2disp_[local_id] = (count++) * m_blockmemsize(1);
        }
    }
    MPI_Win_fence(0, local2disp_window_);

    //................................................
    // post the exposure epoch and start the access one for local2mirrors
    m_assert(mirror_origin_group_ != MPI_GROUP_NULL, "call the InitComm function first!");
    m_assert(mirror_target_group_ != MPI_GROUP_NULL, "call the InitComm function first!");
    MPI_Win_post(mirror_origin_group_, 0, local2disp_window_);
    MPI_Win_start(mirror_target_group_, 0, local2disp_window_);

    // init the list on every active block that matches the level requirements
    for (level_t il = min_level_; il <= max_level_; il++) {
        m_log("init ghosts for level %d", il);
        DoOp_F_<op_t<Ghost*, nullptr_t>, Ghost*, nullptr_t>(CallGhostInitList, grid_, il, nullptr, this);
    }

    // complete the access epoch and wait for the exposure one
    MPI_Win_complete(local2disp_window_);
    MPI_Win_wait(local2disp_window_);

    //................................................
    // local2disp_ = nullptr;
    MPI_Win_free(&local2disp_window_);
    local2disp_window_ = MPI_WIN_NULL;
    //-------------------------------------------------------------------------
    m_end;
}

void Ghost::FreeList_() {
    m_begin;
    //-------------------------------------------------------------------------
    p8est_mesh_t* mesh = grid_->mesh();
    // clear the lists
    for (lid_t ib = 0; ib < n_active_quad_; ib++) {
        // free the blocks
        for (auto biter : (*block_children_[ib])) delete biter;
        for (auto biter : (*ghost_children_[ib])) delete biter;
        for (auto biter : (*block_sibling_[ib])) delete biter;
        for (auto biter : (*ghost_sibling_[ib])) delete biter;
        for (auto biter : (*block_parent_[ib])) delete biter;
        for (auto biter : (*ghost_parent_[ib])) delete biter;
        for (auto piter : (*phys_[ib])) delete piter;
        // purge everything
        delete (block_children_[ib]);
        delete (ghost_children_[ib]);
        delete (block_sibling_[ib]);
        delete (ghost_sibling_[ib]);
        delete (block_parent_[ib]);
        delete (ghost_parent_[ib]);
        delete (phys_[ib]);
    }
    m_free(block_children_);
    m_free(ghost_children_);
    m_free(block_sibling_);
    m_free(ghost_sibling_);
    m_free(block_parent_);
    m_free(ghost_parent_);
    m_free(phys_);
    //-------------------------------------------------------------------------
    m_end;
}

void Ghost::InitComm_(){
    m_begin;
    m_assert(n_active_quad_ >= 0, "the number of active quads must be computed beforehand");
    //-------------------------------------------------------------------------
    // get stupid information
    int            mpi_size = grid_->mpisize();
    MPI_Comm       mpi_comm = grid_->mpicomm();
    p8est_t*       forest   = grid_->forest();
    p8est_ghost_t* ghost    = grid_->ghost();

    //................................................
    // compute the number of admissible local mirrors and store their reference in the array
    n_mirror_to_send_ = 0;
    for (iblock_t im = 0; im < ghost->mirrors.elem_count; im++) {
        qdrt_t* mirror       = p8est_quadrant_array_index(&ghost->mirrors, im);
        level_t   mirror_level = p4est_GetQuadFromMirror(forest, mirror)->level;
        // update the counters if the mirror is admissible
        // we need to have called the function p4est_balance()!!
        if ((min_level_ - 1) <= mirror_level && mirror_level <= (max_level_ + 1)) {
            n_mirror_to_send_++;
        }
    }
    // initialize the Window by allocating the memory space needed for the exchange
    MPI_Info info;
    MPI_Info_create(&info);
    MPI_Info_set(info, "no_locks", "true");
    MPI_Aint win_mem_size = n_mirror_to_send_ * m_blockmemsize(1) * sizeof(real_t);
    MPI_Win_allocate(win_mem_size, sizeof(real_t), info, mpi_comm, &mirrors_, &mirrors_window_);
    MPI_Info_free(&info);

    //................................................
    // get the list of ranks that will generate a call to access my mirrors
    int  n_in_group  = 0;
    int* group_ranks = reinterpret_cast<int*>(m_calloc(mpi_size * sizeof(int)));
    for (int ir = 0; ir < mpi_size; ir++) {
        // for every mirror send to that precise rank
        iblock_t send_first = ghost->mirror_proc_offsets[ir];
        iblock_t send_last  = ghost->mirror_proc_offsets[ir + 1];
        for (iblock_t bid = send_first; bid < send_last; bid++) {
            // access the element and ask for its level
            qdrt_t* mirror       = p8est_quadrant_array_index(&ghost->mirrors, ghost->mirror_proc_mirrors[bid]);
            sid_t   mirror_level = p4est_GetQuadFromMirror(forest, mirror)->level;
            // if the mirror is admissible, register the rank and break
            if ((min_level_ - 1) <= mirror_level && mirror_level <= (max_level_ + 1)) {
                group_ranks[n_in_group++] = ir;
                break;
            }
        }
    }
    // get the RMA mirror group ready
    MPI_Group win_group;
    MPI_Win_get_group(mirrors_window_, &win_group);
    MPI_Group_incl(win_group, n_in_group, group_ranks, &mirror_origin_group_);

    // get the list of ranks that will received a call from me to access their mirrors
    n_in_group = 0;
    for (int ir = 0; ir < mpi_size; ir++) {
        // for every mirror send to that precise rank
        iblock_t recv_first = ghost->proc_offsets[ir];
        iblock_t recv_last  = ghost->proc_offsets[ir + 1];
        for (iblock_t bid = recv_first; bid < recv_last; bid++) {
            // access the element and ask for its level
            p8est_quadrant_t* ghostquad   = p8est_quadrant_array_index(&ghost->ghosts, bid);
            sid_t             ghost_level = ghostquad->level;
            // if the mirror is admissible, register the rank and break
            if ((min_level_ - 1) <= ghost_level && ghost_level <= (max_level_ + 1)) {
                group_ranks[n_in_group++] = ir;
                break;
            }
        }
    }
    MPI_Group_incl(win_group, n_in_group, group_ranks, &mirror_target_group_);
    MPI_Group_free(&win_group);

    //................................................
    // free the allocated memory
    m_free(group_ranks);

    // assert that everybody has created the windows correctly
    MPI_Win_fence(0, mirrors_window_);
    //-------------------------------------------------------------------------
    m_end;
}

void Ghost::FreeComm_() {
    m_begin;
    //-------------------------------------------------------------------------
    // free the group
    MPI_Group_free(&mirror_origin_group_);
    MPI_Group_free(&mirror_target_group_);
    // free the window
    MPI_Win_free(&mirrors_window_);
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief step 1: post the exchange for the coarse and sibling ghost information
 * 
 * @param field 
 * @param ida 
 */
void Ghost::GetGhost_Post(Field* field, const sid_t ida) {
    m_begin;
    m_assert(ida >= 0, "the ida must be >=0!");
    m_assert(grid_->is_mesh_valid(), "the mesh needs to be valid before entering here");
    //-------------------------------------------------------------------------
    // store the current dimension
    ida_ = ida;
    // fill the Window memory with the Mirror information
    LoopOnMirrorBlock_(&Ghost::PushToWindow4Block, field);
    // post the exposure epoch for my own mirrors: I am a target warning that origin group will RMA me
    MPI_Win_post(mirror_origin_group_, 0, mirrors_window_);
    // start the access epoch, to get info from neighbors: I am an origin warning that I will RMA the target group
    MPI_Win_start(mirror_target_group_, 0, mirrors_window_);
    // start what can be done = sibling and parents copy
    for (level_t il = min_level_; il <= max_level_; il++) {
        DoOp_F_<op_t<Ghost*, Field*>, Ghost*, Field*>(CallGetGhost4Block_Post, grid_, il, field, this);
    }
    //-------------------------------------------------------------------------
    m_end;
}

void Ghost::GetGhost_Wait(Field* field, const sid_t ida) {
    m_begin;
    m_assert(ida >= 0, "the ida must be >=0!");
    m_assert(grid_->is_mesh_valid(), "the mesh needs to be valid before entering here");
    //-------------------------------------------------------------------------
    // finish the access epochs for the exposure epoch to be over
    MPI_Win_complete(mirrors_window_);
    // wait for the exposure epoch to be over
    MPI_Win_wait(mirrors_window_);
    m_log("now starting refinement on field %s and dimension %d",field->name().c_str(),ida);
    // we now have all the information needed for the refinement
    for (level_t il = min_level_; il <= max_level_; il++) {
        DoOp_F_<op_t<Ghost*, Field*>, Ghost*, Field*>(CallGetGhost4Block_Wait, grid_, il, field, this);
    }
    //-------------------------------------------------------------------------
    m_end;
}

void Ghost::PutGhost_Post(Field* field, const sid_t ida) {
    m_begin;
    m_assert(ida >= 0, "the ida must be >=0!");
    m_assert(grid_->is_mesh_valid(), "the mesh needs to be valid before entering here");
    //-------------------------------------------------------------------------
    //-------------------------------------------------------------------------
    m_end;
}

void Ghost::PutGhost_Wait(Field* field, const sid_t ida) {
    m_begin;
    m_assert(ida >= 0, "the ida must be >=0!");
    m_assert(grid_->is_mesh_valid(), "the mesh needs to be valid before entering here");
    //-------------------------------------------------------------------------
    //-------------------------------------------------------------------------
    m_end;
}



// /**
//  * @brief given the dimension and the field, copy the mirror blocks that are located on the given level to the send buffers
//  * 
//  * @param field 
//  * @param ida 
//  */
// void Ghost::PushToMirror(Field* field, const sid_t ida) {
//     m_begin;
//     m_assert(grid_->is_mesh_valid(), "the mesh needs to be valid before entering here");
//     //-------------------------------------------------------------------------
//     // // store the current dimension
//     // ida_ = ida;
//     // // loop on the mirrors and copy the values
//     // LoopOnMirrorBlock_(&Ghost::PushToMirror4Block, field);
//     //-------------------------------------------------------------------------
//     m_end;
// }

// /**
//  * @brief starts the send requests for all the mirror block that are located on a given level
//  * 
//  * @warning we do not start the reception requests because we are not sure the buffers are available
//  * 
//  */
// void Ghost::MirrorToGhostSend(Prof* prof) {
//     m_begin;
//     //-------------------------------------------------------------------------
//     // if (n_send_request_ > 0) {
//     //     m_profStart(prof, "ghost_comm_start");
//     //     // youhou
//     //     MPI_Startall(n_send_request_, mirror_send_);
//     //     // finitooo
//     //     m_profStop(prof, "ghost_comm_start");
//     // }
//     //-------------------------------------------------------------------------
//     m_end;
// }

// /**
//  * @brief starts the reception requests, end the send requests and end the reception requests for the blocks that are located on a given level
//  * 
//  * @param prof profiler, might be nullptr
//  */
// void Ghost::MirrorToGhostRecv(Prof* prof) {
//     m_begin;
//     //-------------------------------------------------------------------------
//     // if (n_recv_request_ > 0) {
//     //     m_profStart(prof, "ghost_comm_start");
//     //     //
//     //     MPI_Startall(n_recv_request_, ghost_recv_);
//     //     //
//     //     m_profStop(prof, "ghost_comm_start");
//     // }
//     // if (n_send_request_ > 0) {
//     //     m_profStart(prof, "ghost_comm_wait");
//     //     //
//     //     MPI_Waitall(n_send_request_, mirror_send_, MPI_STATUSES_IGNORE);
//     //     //
//     //     m_profStop(prof, "ghost_comm_wait");
//     // }
//     // if (n_recv_request_ > 0) {
//     //     m_profStart(prof, "ghost_comm_wait");
//     //     //
//     //     MPI_Waitall(n_recv_request_, ghost_recv_, MPI_STATUSES_IGNORE);
//     //     //
//     //     m_profStop(prof, "ghost_comm_wait");
//     // }
//     //-------------------------------------------------------------------------
//     m_end;
// }

// /**
//  * @brief compute the ghost values from the received buffers and the local blocks
//  * 
//  * @param field
//  * @param ida 
//  * @param interp
//  */
// void Ghost::PullFromGhost(Field* field, const sid_t ida) {
//     m_begin;
//     m_assert(grid_->is_mesh_valid(), "the mesh needs to be valid before entering here");
//     //-------------------------------------------------------------------------
//     // // store the current dimension
//     // ida_ = ida;
//     // // interpolate
//     // for (level_t il = min_level_; il <= max_level_; il++) {
//     //     DoOp_F_<op_t<Ghost*, Field*>, Ghost*, Field*>(CallGhostPullFromGhost, grid_, il, field, this);
//     // }
//     //-------------------------------------------------------------------------
//     m_end;
// }

// /**
//  * @brief initalize the MPI_Request communication once, ready to be used
//  * 
//  * @note One p4est mirror structure can be send to multiple ranks
//  * 
//  */
// void Ghost::InitComm_() {
//     m_begin;
//     //-------------------------------------------------------------------------
//     int            mpi_size = grid_->mpisize();
//     MPI_Comm       mpi_comm = grid_->mpicomm();
//     p8est_t*       forest   = grid_->forest();
//     p8est_ghost_t* ghost    = grid_->ghost();

//     // get who send/recv what by level
//     for (int ir = 0; ir < mpi_size; ir++) {
//         // send a mirror
//         iblock_t rank_n_2_send = 0;
//         iblock_t send_first    = ghost->mirror_proc_offsets[ir];
//         iblock_t send_last     = ghost->mirror_proc_offsets[ir + 1];
//         for (iblock_t bid = send_first; bid < send_last; bid++) {
//             // access the element and ask for its level
//             qdrt_t* mirror       = p8est_quadrant_array_index(&ghost->mirrors, ghost->mirror_proc_mirrors[bid]);
//             sid_t   mirror_level = p4est_GetQuadFromMirror(forest, mirror)->level;
//             // update the counters if the mirror is admissible
//             // we need to have called the function p4est_balance()!!
//             if ((min_level_ - 1) <= mirror_level && mirror_level <= (max_level_ + 1)) {
//                 rank_n_2_send++;
//             }
//         }
//         // if we added some blocks, create a new request and increment the global counter
//         n_send_request_ += (rank_n_2_send > 0) ? 1 : 0;
//         n_mirror_to_send_ += rank_n_2_send;

//         // recv a ghost
//         iblock_t rank_n_2_recv = 0;
//         iblock_t recv_first    = ghost->proc_offsets[ir];
//         iblock_t recv_last     = ghost->proc_offsets[ir + 1];
//         for (iblock_t bid = recv_first; bid < recv_last; bid++) {
//             // get the ghost level
//             level_t ghost_level = p8est_quadrant_array_index(&ghost->ghosts, bid)->level;
//             // it is not ensured that the ghost quad has a level, so run a check
//             if ((min_level_ - 1) <= ghost_level && ghost_level <= (max_level_ + 1)) {
//                 rank_n_2_recv++;
//             }
//         }
//         // if we added some blocks, create a new request and increment the global counter
//         n_recv_request_ += (rank_n_2_recv > 0) ? 1 : 0;
//         n_ghost_to_recv_ += rank_n_2_recv;
//     }

//     // allocate the request and local arrays for each level
//     local_to_mirrors = reinterpret_cast<iblock_t*>(m_calloc(n_mirror_to_send_ * sizeof(iblock_t*)));
//     ghost_to_local_  = reinterpret_cast<iblock_t*>(m_calloc(ghost->ghosts.elem_count * sizeof(iblock_t*)));
//     mirror_send_     = reinterpret_cast<MPI_Request*>(m_calloc(n_send_request_ * sizeof(MPI_Request)));
//     ghost_recv_      = reinterpret_cast<MPI_Request*>(m_calloc(n_recv_request_ * sizeof(MPI_Request)));

//     // allocate the mirror and ghost data arrays to the max among each level
//     mirrors_ = reinterpret_cast<real_t*>(m_calloc(sizeof(real_t) * M_NGHOST * n_mirror_to_send_));
//     ghosts_  = reinterpret_cast<real_t*>(m_calloc(sizeof(real_t) * M_NGHOST * n_ghost_to_recv_));

//     // allocate the MPI requests
//     lid_t    send_request_offset = 0;
//     lid_t    recv_request_offset = 0;
//     iblock_t mirror_offset       = 0;
//     iblock_t ghost_offset        = 0;
//     for (int ir = 0; ir < mpi_size; ir++) {
//         // send a mirror
//         iblock_t rank_n_2_send = 0;
//         iblock_t send_first    = ghost->mirror_proc_offsets[ir];
//         iblock_t send_last     = ghost->mirror_proc_offsets[ir + 1];
//         for (iblock_t bid = send_first; bid < send_last; bid++) {
//             // access the element and ask for its level
//             qdrt_t* mirror       = p8est_quadrant_array_index(&ghost->mirrors, ghost->mirror_proc_mirrors[bid]);
//             sid_t   mirror_level = p4est_GetQuadFromMirror(forest, mirror)->level;
//             if ((min_level_ - 1) <= mirror_level && mirror_level <= (max_level_ + 1)) {
//                 iblock_t add_id          = mirror_offset + rank_n_2_send;
//                 local_to_mirrors[add_id] = ghost->mirror_proc_mirrors[bid];
//                 rank_n_2_send++;
//             }
//         }
//         if (rank_n_2_send > 0) {
//             // move the sendbuf to howmany have already been send * the size of one block
//             real_p send_buf  = mirrors_ + mirror_offset * M_NGHOST;
//             size_t send_size = rank_n_2_send * M_NGHOST;
//             MPI_Send_init(send_buf, send_size, M_MPI_REAL, ir, P4EST_COMM_GHOST_EXCHANGE, mpi_comm, mirror_send_ + send_request_offset);
//             // update the offsets
//             mirror_offset += rank_n_2_send;
//             send_request_offset++;
//         }
//         // recv a ghost
//         iblock_t rank_n_2_recv = 0;
//         iblock_t recv_first    = ghost->proc_offsets[ir];
//         iblock_t recv_last     = ghost->proc_offsets[ir + 1];
//         for (iblock_t bid = recv_first; bid < recv_last; bid++) {
//             // get the ghost level
//             level_t ghost_level = p8est_quadrant_array_index(&ghost->ghosts, bid)->level;
//             // it is not ensured that the ghost quad has a level, so run a check
//             if ((min_level_ - 1) <= ghost_level && ghost_level <= (max_level_ + 1)) {
//                 iblock_t add_id      = ghost_offset + rank_n_2_recv;
//                 ghost_to_local_[bid] = add_id;
//                 rank_n_2_recv++;
//             }
//         }
//         if (rank_n_2_recv > 0) {
//             // move the sendbuf to howmany have already been send * the size of one block
//             real_p recv_buf  = ghosts_ + ghost_offset * M_NGHOST;
//             size_t recv_size = rank_n_2_recv * M_NGHOST;
//             MPI_Recv_init(recv_buf, recv_size, M_MPI_REAL, ir, P4EST_COMM_GHOST_EXCHANGE, mpi_comm, ghost_recv_ + recv_request_offset);
//             // update the offsets
//             ghost_offset += rank_n_2_recv;
//             recv_request_offset++;
//         }
//     }
//     m_assert(ghost_offset == n_ghost_to_recv_, "the two numbers have to match");
//     m_assert(mirror_offset == n_mirror_to_send_, "the two numbers have to match");
//     m_assert(send_request_offset == n_send_request_, "the two numbers have to match");
//     m_assert(recv_request_offset == n_recv_request_, "the two numbers have to match");

//     //-------------------------------------------------------------------------
//     m_end;
// }

/**
 * @brief setup the ghosting lists, for the considered block
 * 
 * @warning this function is called within a omp parellel region.
 * The p4est arrays are NOT thread-safe, while the lists have been assumed to be NOT threadsafe as well.
 * 
 * @note the p4est array allocation is thread-safe, while the tracking of the memory allocated is NOT.
 * 
 * @param qid the quadrant ID
 * @param block the associated block
 */
void Ghost::InitList4Block(const qid_t* qid, GridBlock* block) {
    //-------------------------------------------------------------------------
    // get the current lists
    ListGBLocal*  bchildren = block_children_[qid->cid];
    ListGBLocal*  bsibling  = block_sibling_[qid->cid];
    ListGBLocal*  bparent   = block_parent_[qid->cid];
    ListGBMirror* gchildren = ghost_children_[qid->cid];
    ListGBMirror* gsibling  = ghost_sibling_[qid->cid];
    ListGBMirror* gparent   = ghost_parent_[qid->cid];
    listGBPhysic* phys      = phys_[qid->cid];

    //................................................
    // allocate the ghost pointer
    block->AllocatePtrGhost(CoarseMemSize(interp_));

    //................................................
    // temporary sc array used to get the ghosts
    sc_array_t* ngh_quad;  // points to the quad
    sc_array_t* ngh_qid;   // give the ID of the quad or ghost
    sc_array_t* ngh_enc;   // get the status

#pragma omp critical
    {
        ngh_quad = sc_array_new(sizeof(qdrt_t*));
        ngh_enc  = sc_array_new(sizeof(int));
        ngh_qid  = sc_array_new(sizeof(int));
    }

    p8est_t*       forest = grid_->forest();
    p8est_mesh_t*  mesh   = grid_->mesh();
    p8est_ghost_t* ghost  = grid_->ghost();

    for (sid_t ibidule = 0; ibidule < M_NNEIGHBOR; ibidule++) {
        //................................................
#pragma omp critical
        {
            // get the neighboring quadrant
            sc_array_reset(ngh_quad);
            sc_array_reset(ngh_enc);
            sc_array_reset(ngh_qid);
            // m_log("I have a grid with %d blocks and %d ghosts. Looking for neighbor %d of block %d = %d, %d",mesh->local_num_quadrants,ghost->ghosts.elem_count,ibidule,qid->cid,qid->qid,qid->tid);
            p8est_mesh_get_neighbors(forest, ghost, mesh, qid->cid, ibidule, ngh_quad, ngh_enc, ngh_qid);
        }
        const iblock_t nghosts = ngh_enc->elem_count;

        //................................................
        // get the number of ghost:
        sid_t nghost_front[3];
        sid_t nghost_back[3];
        for (sid_t id = 0; id < 3; id++) {
            // set the number of ghost to compute
            nghost_front[id] = interp_->nghost_front();
            nghost_back[id]  = interp_->nghost_back();
        }

        //................................................
        // no ghosts? then is a physical BC
        if (nghosts == 0) {
            // we only apply the physics to entire faces
            if (ibidule < 6) {
                PhysBlock* pb = new PhysBlock(ibidule, block, nghost_front, nghost_back);
#pragma omp critical
                phys->push_back(pb);
            }
            // else, the edges and corners will be filled through the face
        }

        //................................................
        // this is a real block or a ghost
        for (iblock_t nid = 0; nid < nghosts; nid++) {
            const int  status  = *(ngh_enc->array + nid * sizeof(int));
            const bool isghost = (status < 0);
            qdrt_t*    nghq    = p4est_GetElement<qdrt_t*>(ngh_quad, nid);

            // get the sign, i.e. the normal to the face, the edge of the corner we consider
            real_t sign[3];
            GhostGetSign(ibidule, sign);

            //................................................
            // get the position of the neighbor, as seen by me!!! may be different than the actual one if there is a periodic bc
            real_t ngh_pos[3];
            if (!isghost) {
                // cannot use the p8est function because the which_tree is not assigned, so we retrieve the position through the block
                GridBlock* ngh_block = *(reinterpret_cast<GridBlock**>(nghq->p.user_data));
                ngh_pos[0]           = ngh_block->xyz(0);
                ngh_pos[1]           = ngh_block->xyz(1);
                ngh_pos[2]           = ngh_block->xyz(2);
            } else {
                p8est_qcoord_to_vertex(grid_->connect(), nghq->p.piggy3.which_tree, nghq->x, nghq->y, nghq->z, ngh_pos);
            }
            // fix the shift in coordinates needed IF the domain is periodic
            for (sid_t id = 0; id < 3; id++) {
                // if we are periodic, we overwrite the position in the direction of the normal !!ONLY!!
                // since it is my neighbor in this normal direction, I am 100% sure that it's origin corresponds to the end of my block
                const real_t to_replace = sign[id] * sign[id] * grid_->domain_periodic(id);  // is (+-1)^2 = +1 if we need to replace it, 0.0 otherwize
                // get the expected position
                const real_t expected_pos = block->xyz(id) + (sign[id] > 0.5) * m_quad_len(block->level()) - (sign[id] < -0.5) * m_quad_len(nghq->level);
                // we override the position if a replacement is needed only
                ngh_pos[id] = to_replace * expected_pos + (1.0 - to_replace) * ngh_pos[id];
            }

            //................................................
            // create the new block and push back
            if (!isghost) {
                // create the new local block
                GBLocal* gb = new GBLocal(block, nghq->level, ngh_pos, nghost_front, nghost_back, true);

                // associate the corresponding neighboring block
                GridBlock* ngh_block = *(reinterpret_cast<GridBlock**>(nghq->p.user_data));
                gb->data_src(ngh_block);

                // register the gb in a list
                if (gb->dlvl() == 0) {
#pragma omp critical
                    bsibling->push_back(gb);
                } else if (gb->dlvl() == -1) {
#pragma omp critical
                    bparent->push_back(gb);
                    // get my contribution to the ghost of my parent neighbor
                    GBLocal* invert_gb = new GBLocal(block, nghq->level, ngh_pos, nghost_front, nghost_back, false);
                    invert_gb->data_src(ngh_block);
#pragma omp critical
                    bchildren->push_back(invert_gb);
                } else {
                    m_assert(gb->dlvl() == 1, "The delta level is not correct: %d", gb->dlvl());
                }
            } else {
                // get the local number in the remote rank and the remote rank
                rank_t ngh_local_id = nghq->p.piggy3.local_num;
                rank_t ngh_rank     = p4est_GetOwnerFromGhost(forest, nghq);
                m_assert(ngh_rank > -1, "p4est unable to recover the rank... baaaad news");

                // create the new mirror block
                GBMirror* gb = new GBMirror(block, nghq->level, ngh_pos, nghost_front, nghost_back, true, ngh_rank);
                // ask the displacement (will be available later, when completing the call
                MPI_Get(gb->data_src_ptr(), 1, MPI_AINT, ngh_rank, ngh_local_id, 1, MPI_AINT, local2disp_window_);

                // register the ghost block in a list
                if (gb->dlvl() == 0) {
#pragma omp critical
                    gsibling->push_back(gb);
                } else if (gb->dlvl() == -1) {
#pragma omp critical
                    gparent->push_back(gb);
                    // get my contribution to the ghost of my parent neighbor
                    GBLocal* invert_gb = new GBLocal(block, nghq->level, ngh_pos, nghost_front, nghost_back, false, ngh_rank);
                    MPI_Get(invert_gb->data_src(), 1, MPI_AINT, ngh_rank, ngh_local_id, 1, MPI_AINT, local2disp_window_);
#pragma omp critical
                    bchildren->push_back(invert_gb);
                } else {
                    m_assert(gb->dlvl() == 1, "The delta level is not correct: %d", gb->dlvl());
                }
            }
        }
    }
#pragma omp critical
    {
        sc_array_destroy(ngh_quad);
        sc_array_destroy(ngh_enc);
        sc_array_destroy(ngh_qid);
    }
    //-------------------------------------------------------------------------
}

// /**
//  * @brief copy the block information into the mirror buffer for a given field and dimension
//  * 
//  * @param qid the quadrant ID
//  * @param block th block that is an actual mirror
//  * @param fid the field ID
//  */
// void Ghost::PushToMirror4Block(const qid_t* qid, GridBlock* block, Field* fid) {
//     m_assert(ida_ >= 0, "the current working dimension has to be correct");
//     m_assert(ida_ < fid->lda(), "the current working dimension has to be correct");
//     //-------------------------------------------------------------------------
//     real_p mirror = mirrors_ + qid->cid * M_NGHOST;
//     real_p data   = block->data(fid, ida_);

//     for (int i2 = 0; i2 < M_N; i2++) {
//         for (int i1 = 0; i1 < M_N; i1++) {
//             for (int i0 = 0; i0 < M_N; i0++) {
//                 mirror[m_sidx(i0, i1, i2, 0, M_N)] = data[m_idx(i0, i1, i2)];
//             }
//         }
//     }
//     //-------------------------------------------------------------------------
// }

void Ghost::PushToWindow4Block(const qid_t* qid, GridBlock* block, Field* fid) {
    m_assert(ida_ >= 0, "the current working dimension has to be correct");
    m_assert(ida_ < fid->lda(), "the current working dimension has to be correct");
    //-------------------------------------------------------------------------
    real_p mirror = mirrors_ + qid->cid * m_blockmemsize(1);
    real_p data   = block->pointer(fid, ida_);
    // data should be aligned
    m_assume_aligned(mirror);
    m_assume_aligned(data);
    // copy the data
    memcpy(data, mirror, m_blockmemsize(1) * sizeof(real_t));
    //-------------------------------------------------------------------------
}

void Ghost::GetGhost4Block_Post(const qid_t* qid, GridBlock* block, Field* fid) {
    //-------------------------------------------------------------------------
    // get the working array given the thread
    real_p tmp = block->ptr_ghost();
    // determine if we have to use the coarse representation
    const bool do_coarse = (block_parent_[qid->cid]->size() + ghost_parent_[qid->cid]->size()) > 0;
    // if so, reset the coarse info
    if (do_coarse) {
        memset(tmp, 0, CoarseMemSize(interp_));
    }

    if (qid->cid == 16) {
        lid_t coarse_start[3];
        lid_t coarse_end[3];

        //set the coarse block to its whole domain
        for (int id = 0; id < 3; id++) {
            coarse_start[id] = 0;
            coarse_end[id]   = M_HN;
        }
        SubBlock* coarse_block = new SubBlock(CoarseNGhostFront(interp_), CoarseStride(interp_), coarse_start, coarse_end);
        for (int i2 = 7; i2 < 8; i2++) {
            for (int i1 = 0; i1 < 9; i1++) {
                printf("%d %d: ", i1, i2);
                for (int i0 = 0; i0 < M_HN; i0++) {
                    printf("%f ", *(tmp + m_zeroidx(0, coarse_block) + m_sidx(i0, i1, i2, 0, coarse_block->stride())));
                }
                printf("\n");
            }
        }
    }

    //................................................
    // start with the ghost blocks
    Compute4Block_GetRma2Myself_(ghost_sibling_[qid->cid], fid, block, block->data(fid, ida_));
    if (do_coarse) {
        Compute4Block_GetRma2Coarse_(ghost_sibling_[qid->cid], fid, block, tmp);
        Compute4Block_GetRma2Coarse_(ghost_parent_[qid->cid], fid, block, tmp);
    }

    // and now the local ones!
    Compute4Block_Copy2Myself_(block_sibling_[qid->cid], fid, block, block->data(fid, ida_));
    if (do_coarse) {
        m_log("-----------------");
        m_log("QID = %d", qid->cid);
        m_log("doing copy to coarse");
        Compute4Block_Copy2Coarse_(block_sibling_[qid->cid], fid, block, tmp);
        Compute4Block_Copy2Coarse_(block_parent_[qid->cid], fid, block, tmp);
    }
    if (qid->cid == 16) {
        m_log("check the coarse after copying the parents and the siblings");
        lid_t coarse_start[3];
        lid_t coarse_end[3];

        //set the coarse block to its whole domain
        for (int id = 0; id < 3; id++) {
            coarse_start[id] = 0;
            coarse_end[id]   = M_HN;
        }
        SubBlock* coarse_block = new SubBlock(CoarseNGhostFront(interp_), CoarseStride(interp_), coarse_start, coarse_end);
        for (int i2 = -1; i2 < 9; i2++) {
            for (int i1 = 4; i1 < 5; i1++) {
                printf("x %d %d: ", i1, i2);
                for (int i0 = -1; i0 < 1; i0++) {
                    printf("%f ", *(tmp + m_zeroidx(0, coarse_block) + m_sidx(i0, i1, i2, 0, coarse_block->stride())));
                }
                printf("\n");
            }
        }
    }

    //-------------------------------------------------------------------------
}

void Ghost::GetGhost4Block_Wait(const qid_t* qid, GridBlock* block, Field* fid) {
    //-------------------------------------------------------------------------
    // get the working array given the thread
    real_p tmp = block->ptr_ghost();
    // determine if we have to use the coarse representation
    const bool do_coarse = (block_parent_[qid->cid]->size() + ghost_parent_[qid->cid]->size()) > 0;

    if (qid->cid == 16) {
        m_log("check the coarse before copying myself");
        lid_t coarse_start[3];
        lid_t coarse_end[3];

        //set the coarse block to its whole domain
        for (int id = 0; id < 3; id++) {
            coarse_start[id] = 0;
            coarse_end[id]   = M_HN;
        }
        SubBlock* coarse_block = new SubBlock(CoarseNGhostFront(interp_), CoarseStride(interp_), coarse_start, coarse_end);
        for (int i2 = -1; i2 < 9; i2++) {
            for (int i1 = 4; i1 < 5; i1++) {
                printf("x %d %d: ", i1, i2);
                for (int i0 = -1; i0 < 1; i0++) {
                    printf("%f ", *(tmp + m_zeroidx(0, coarse_block) + m_sidx(i0, i1, i2, 0, coarse_block->stride())));
                }
                printf("\n");
            }
        }
    }

    // now that the copy are over, just need to add myself and the physics
    if (do_coarse) {
        Compute4Block_Myself2Coarse_(qid, block, fid, tmp);

        if (qid->cid == 16) {
            m_log("check the coarse after copying myself");
            lid_t coarse_start[3];
            lid_t coarse_end[3];

            //set the coarse block to its whole domain
            for (int id = 0; id < 3; id++) {
                coarse_start[id] = 0;
                coarse_end[id]   = M_HN;
            }
            SubBlock* coarse_block = new SubBlock(CoarseNGhostFront(interp_), CoarseStride(interp_), coarse_start, coarse_end);
            for (int i2 = -1; i2 < 9; i2++) {
                for (int i1 = 4; i1 < 5; i1++) {
                    printf("x %d %d: ",i1,i2);
                    for (int i0 = -1; i0 < 1; i0++) {
                        printf("%f ", *(tmp + m_zeroidx(0, coarse_block) + m_sidx(i0, i1, i2, 0, coarse_block->stride())));
                    }
                    printf("\n");
                }
            }
        }
        // compute the fine ghosts
        m_log("-----------------");
        m_log("QID = %d", qid->cid);
        Compute4Block_Refine_(block_parent_[qid->cid], tmp, block->data(fid, ida_));
        Compute4Block_Refine_(ghost_parent_[qid->cid], tmp, block->data(fid, ida_));
    } else {
        m_log("nothing to do, leaving");
    }

    //-------------------------------------------------------------------------
}

void Ghost::PutGhost4Block_Post(const qid_t *qid, GridBlock *block, Field *fid){

}
void Ghost::PutGhost4Block_Wait(const qid_t *qid, GridBlock *block, Field *fid){

}

// /**
//  * @brief given an up-to-date @ref ghost_ arrays,
//  * compute the ghosts points for the current block (field @ref fid, dimension @ref ida_)
//  * 
//  * @warning this function runs inside an OpenMP parallel region
//  * 
//  * @note a lot of dependencies appear if the block has a coarser neighbor. To solve it, we build
//  * step by step a coarse representation of the current block, which is then used to refine to the actual block.
//  * 
//  * @param qid the quadrant id, see @ref qid_t
//  * @param cur_block the current grid block, see @ref GridBlock
//  * @param fid the current field, see @ref Field
//  */
// void Ghost::PullFromGhost4Block(const qid_t* qid, GridBlock* cur_block, Field* fid) {
//     //-------------------------------------------------------------------------
//     // get a subblock describing the ghost memory
//     lid_t     ghost_start[3] = {0, 0, 0};
//     lid_t     ghost_end[3]   = {M_N, M_N, M_N};
//     SubBlock* ghost_subblock = new SubBlock(0, M_N, ghost_start, ghost_end);

//     lid_t     coarse_start[3] = {0, 0, 0};
//     lid_t     coarse_end[3]   = {M_HN, M_HN, M_HN};
//     SubBlock* coarse_subblock = new SubBlock(CoarseNGhostFront(interp_), CoarseStride(interp_), coarse_start, coarse_end);

//     // get the working array given the thread
//     const int ithread = omp_get_thread_num();
//     real_p    tmp     = coarse_tmp_[ithread];

//     // determine if we have to use the coarse representation
//     const bool do_coarse = (block_parent_[qid->cid]->size() + ghost_parent_[qid->cid]->size()) > 0;
//     // if so, reset the coarse info
//     if (do_coarse) {
//         memset(tmp, 0, CoarseStride(interp_) * CoarseStride(interp_) * CoarseStride(interp_) * sizeof(real_t));
//     }

//     // first, treat the siblings -> will copy to the ghost location + copy to the temp
//     PullFromGhost4Block_Sibling_(qid, cur_block, fid, do_coarse, ghost_subblock, coarse_subblock, tmp);

//     // do a coarse version of myself and complete with some physics if needed
//     if (do_coarse) {
//         PullFromGhost4Block_Myself_(qid, cur_block, fid, coarse_subblock, tmp);
//         // reset the coarse sublock to the full position
//         for (int id = 0; id < 3; id++) {
//             coarse_start[id] = -CoarseNGhostFront(interp_);
//             coarse_end[id]   = CoarseStride(interp_) - CoarseNGhostFront(interp_);
//         }
//         coarse_subblock->Reset(CoarseNGhostFront(interp_), CoarseStride(interp_), coarse_start, coarse_end);
//     }
//     // then, treat the parents -> copy to the temp
//     PullFromGhost4Block_FromParent_(qid, cur_block, fid, do_coarse, ghost_subblock, coarse_subblock, tmp);

//     // //-------------------------------------------------------------------------
//     // // do the blocks, on my level or finer
//     // for (auto biter = block_sibling_[qid->cid]->begin(); biter != block_sibling_[qid->cid]->end(); biter++) {
//     //     GhostBlock* gblock    = (*biter);
//     //     GridBlock*  ngh_block = gblock->block_src();
//     //     // memory details
//     //     MemLayout* block_src = ngh_block;
//     //     real_p     data_src  = ngh_block->data(fid, ida_);
//     //     MemLayout* block_trg = gblock;
//     //     real_p     data_trg  = cur_block->data(fid, ida_);
//     //     // interpolate
//     //     interp_->Interpolate(gblock->dlvl(), gblock->shift(), block_src, data_src, block_trg, data_trg);

//     //     // we need to interpolate on the coarse version of myself as well
//     //     if (do_coarse) {
//     //         // set the coarse block to the correct position
//     //         for (int id = 0; id < 3; id++) {
//     //             coarse_start[id] = CoarseFromBlock(gblock->start(id), cgs(interp_));
//     //             coarse_end[id]   = CoarseFromBlock(gblock->end(id), cgs(interp_));
//     //         }
//     //         coarse_subblock->Reset(cgs(interp_), CoarseStride(interp_), coarse_start, coarse_end);
//     //         // memory details
//     //         MemLayout* block_src = ngh_block;
//     //         real_p     data_src  = ngh_block->data(fid, ida_);
//     //         MemLayout* block_trg = coarse_subblock;
//     //         real_p     data_trg  = tmp + m_zeroidx(0, coarse_subblock);
//     //         // interpolate, the level is 1 coarser and the shift is unchanged
//     //         interp_->Interpolate(gblock->dlvl() + 1, gblock->shift(), block_src, data_src, block_trg, data_trg);
//     //     }
//     // }

//     // // do the ghosts, on my level or finer
//     // for (auto biter = ghost_sibling_[qid->cid]->begin(); biter != ghost_sibling_[qid->cid]->end(); biter++) {
//     //     GhostBlock* gblock = (*biter);
//     //     // memory details
//     //     MemLayout* block_src = ghost_subblock;
//     //     real_p     data_src  = gblock->data_src();
//     //     MemLayout* block_trg = gblock;
//     //     real_p     data_trg  = cur_block->data(fid, ida_);
//     //     // interpolate
//     //     interp_->Interpolate(gblock->dlvl(), gblock->shift(), block_src, data_src, block_trg, data_trg);

//     //     // we need to interpolate on the coarse version of myself as well
//     //     if (do_coarse) {
//     //         // set the coarse block to the correct position
//     //         for (int id = 0; id < 3; id++) {
//     //             coarse_start[id] = CoarseFromBlock(gblock->start(id), cgs(interp_));
//     //             coarse_end[id]   = CoarseFromBlock(gblock->end(id), cgs(interp_));
//     //         }
//     //         coarse_subblock->Reset(cgs(interp_), CoarseStride(interp_), coarse_start, coarse_end);
//     //         // memory details
//     //         MemLayout* block_src = ghost_subblock;
//     //         real_p     data_src  = gblock->data_src();
//     //         MemLayout* block_trg = coarse_subblock;
//     //         real_p     data_trg  = tmp + m_zeroidx(0, coarse_subblock);
//     //         // interpolate, the level is 1 coarser and the shift is unchanged
//     //         interp_->Interpolate(gblock->dlvl() + 1, gblock->shift(), block_src, data_src, block_trg, data_trg);
//     //     }
//     // }

//     // //-------------------------------------------------------------------------
//     // // copy the coarse blocks to the coarse representation
//     // for (auto biter = block_parent_[qid->cid]->begin(); biter != block_parent_[qid->cid]->end(); biter++) {
//     //     GhostBlock* gblock    = (*biter);
//     //     GridBlock*  ngh_block = gblock->block_src();
//     //     // setup the coarse sublock to the position
//     //     for (int id = 0; id < 3; id++) {
//     //         coarse_start[id] = CoarseFromBlock(gblock->start(id), cgs(interp_));
//     //         coarse_end[id]   = CoarseFromBlock(gblock->end(id), cgs(interp_));
//     //     }
//     //     coarse_subblock->Reset(cgs(interp_), CoarseStride(interp_), coarse_start, coarse_end);
//     //     // memory details
//     //     MemLayout* block_src = ngh_block;
//     //     real_p     data_src  = ngh_block->data(fid, ida_);
//     //     MemLayout* block_trg = coarse_subblock;
//     //     real_p     data_trg  = tmp + m_zeroidx(0, coarse_subblock);
//     //     // interpolate, the level is 1 coarser and the shift is unchanged
//     //     m_assert(gblock->dlvl() + 1 == 0, "the gap in level has to be 0");
//     //     interp_->Interpolate(gblock->dlvl() + 1, gblock->shift(), block_src, data_src, block_trg, data_trg);
//     // }
//     // // copy the ghost into the coarse representation
//     // for (auto biter = ghost_parent_[qid->cid]->begin(); biter != ghost_parent_[qid->cid]->end(); biter++) {
//     //     GhostBlock* gblock = (*biter);
//     //     // update the coarse subblock
//     //     for (int id = 0; id < 3; id++) {
//     //         coarse_start[id] = CoarseFromBlock(gblock->start(id), cgs(interp_));
//     //         coarse_end[id]   = CoarseFromBlock(gblock->end(id), cgs(interp_));
//     //     }
//     //     coarse_subblock->Reset(cgs(interp_), CoarseStride(interp_), coarse_start, coarse_end);
//     //     // memory details
//     //     MemLayout* block_src = ghost_subblock;
//     //     real_p     data_src  = gblock->data_src();
//     //     MemLayout* block_trg = coarse_subblock;
//     //     real_p     data_trg  = tmp + m_zeroidx(0, coarse_subblock);
//     //     // interpolate, the level is 1 coarser and the shift is unchanged
//     //     m_assert(gblock->dlvl() + 1 == 0, "the gap in level has to be 0");
//     //     interp_->Interpolate(gblock->dlvl() + 1, gblock->shift(), block_src, data_src, block_trg, data_trg);
//     // }

//     // m_log("dbg: tree %d, quad %d: doing coarse ? %d", qid->tid, qid->qid,do_coarse);
//     //-------------------------------------------------------------------------
//     // interpolate the ghost representation on myself
    
//     // need to put the ghost to my parent's place
//     PullFromGhost4Block_ToParent_(qid, cur_block, fid, ghost_subblock, coarse_subblock, tmp);

//     //-------------------------------------------------------------------------
//     // finally do some physics
//     for (auto piter = phys_[qid->cid]->begin(); piter != phys_[qid->cid]->end(); piter++) {
//         PhysBlock* gblock = (*piter);
//         // get the direction and the corresponding bctype
//         bctype_t bctype = fid->bctype(ida_, gblock->iface());
//         if (bctype == M_BC_EVEN) {
//             EvenBoundary_4 bc = EvenBoundary_4();
//             bc(gblock->iface(), face_start[gblock->iface()], cur_block->hgrid(), 0.0, gblock, cur_block->data(fid, ida_));
//         } else if (bctype == M_BC_ODD) {
//             OddBoundary_4 bc = OddBoundary_4();
//             bc(gblock->iface(), face_start[gblock->iface()], cur_block->hgrid(), 0.0, gblock, cur_block->data(fid, ida_));
//         } else if (bctype == M_BC_EXTRAP_3) {
//             ExtrapBoundary_3 bc = ExtrapBoundary_3();
//             bc(gblock->iface(), face_start[gblock->iface()], cur_block->hgrid(), 0.0, gblock, cur_block->data(fid, ida_));
//         } else if (bctype == M_BC_EXTRAP_4) {
//             ExtrapBoundary_4 bc = ExtrapBoundary_4();
//             bc(gblock->iface(), face_start[gblock->iface()], cur_block->hgrid(), 0.0, gblock, cur_block->data(fid, ida_));
//         } else if (bctype == M_BC_EXTRAP_5) {
//             ExtrapBoundary_5 bc = ExtrapBoundary_5();
//             bc(gblock->iface(), face_start[gblock->iface()], cur_block->hgrid(), 0.0, gblock, cur_block->data(fid, ida_));
//         } else if (bctype == M_BC_ZERO) {
//             ZeroBoundary bc = ZeroBoundary();
//             bc(gblock->iface(), face_start[gblock->iface()], cur_block->hgrid(), 0.0, gblock, cur_block->data(fid, ida_));
//         } else {
//             m_assert(false, "this type of BC is not implemented yet or not valid %d", bctype);
//         }
//     }
//     delete (ghost_subblock);
//     delete (coarse_subblock);

//     // m_log("dbg: tree %d, quad %d: finiiiish", qid->tid, qid->qid);
//     //-------------------------------------------------------------------------
// }

// void Ghost::PullFromGhost4Block_Children(const qid_t *qid, GridBlock *cur_block, Field *fid,
//                                               const bool do_coarse, SubBlock *ghost_block, SubBlock *coarse_block, real_t *coarse_mem){
    // //-------------------------------------------------------------------------
    // lid_t coarse_start[3];
    // lid_t coarse_end[3];

    // // do the blocks on my level
    // for (auto biter = block_children_[qid->cid]->begin(); biter != block_children_[qid->cid]->end(); biter++) {
    //     GhostBlock* gblock    = (*biter);
    //     GridBlock*  ngh_block = gblock->block_src();
    //     // memory details
    //     MemLayout* block_src = ngh_block;
    //     real_p     data_src  = ngh_block->data(fid, ida_);
    //     MemLayout* block_trg = gblock;
    //     real_p     data_trg  = cur_block->data(fid, ida_);
    //     // interpolate
    //     m_assert(gblock->dlvl() == +1,"we are treating children, the delta level must be +1");
    //     interp_->Interpolate(gblock->dlvl(), gblock->shift(), block_src, data_src, block_trg, data_trg);

    //     // we need to copy the data to the coarse representation
    //     if (do_coarse) {
    //         // set the coarse block to the correct position
    //         for (int id = 0; id < 3; id++) {
    //             coarse_start[id] = CoarseFromBlock(gblock->start(id), cgs(interp_));
    //             coarse_end[id]   = CoarseFromBlock(gblock->end(id), cgs(interp_));
    //         }
    //         coarse_block->Reset(cgs(interp_), CoarseStride(interp_), coarse_start, coarse_end);
    //         // memory details
    //         MemLayout* block_src = ngh_block;
    //         real_p     data_src  = ngh_block->data(fid, ida_);
    //         MemLayout* block_trg = coarse_block;
    //         real_p     data_trg  = coarse_mem + m_zeroidx(0, coarse_block);
    //         // interpolate, the level is 1 coarser and the shift is unchanged
    //         m_assert(gblock->dlvl()==1,"the difference of level MUST be 1");
    //         interp_->Copy(gblock->dlvl() + 1, gblock->shift(), block_src, data_src, block_trg, data_trg);
    //     }
    // }

    // // do the ghosts, on my level
    // for (auto biter = ghost_children_[qid->cid]->begin(); biter != ghost_children_[qid->cid]->end(); biter++) {
    //     GhostBlock* gblock = (*biter);
    //     // memory details
    //     MemLayout* block_src = ghost_block;
    //     real_p     data_src  = gblock->data_src();
    //     MemLayout* block_trg = gblock;
    //     real_p     data_trg  = cur_block->data(fid, ida_);
    //     // interpolate
    //     m_assert(gblock->dlvl() == 0,"we are treating sibling, the delta level must be 0");
    //     interp_->Interpolate(gblock->dlvl(), gblock->shift(), block_src, data_src, block_trg, data_trg);

//     // we need to interpolate on the coarse version of myself as well
//     if (do_coarse) {
//         // set the coarse block to the correct position
//         for (int id = 0; id < 3; id++) {
//             coarse_start[id] = CoarseFromBlock(gblock->start(id), cgs(interp_));
//             coarse_end[id]   = CoarseFromBlock(gblock->end(id), cgs(interp_));
//         }
//         coarse_block->Reset(cgs(interp_), CoarseStride(interp_), coarse_start, coarse_end);
//         // memory details
//         MemLayout* block_src = ghost_block;
//         real_p     data_src  = gblock->data_src();
//         MemLayout* block_trg = coarse_block;
//         real_p     data_trg  = coarse_mem + m_zeroidx(0, coarse_block);
//         // interpolate, the level is 1 coarser and the shift is unchanged
//         m_assert(gblock->dlvl()==1,"the difference of level MUST be 1");
//         interp_->Interpolate(gblock->dlvl() + 1, gblock->shift(), block_src, data_src, block_trg, data_trg);
//     }
// }
// //-------------------------------------------------------------------------
// }

void Ghost::Compute4Block_Copy2Myself_(const ListGBLocal* ghost_list, Field* fid, GridBlock* block_trg, real_t* data_trg) {
    //-------------------------------------------------------------------------
    // the source block is always the same
    const lid_t src_start[3] = {0, 0, 0};
    const lid_t src_end[3]   = {M_N, M_N, M_N};
    MemLayout*  block_src    = new SubBlock(M_GS, M_STRIDE, src_start, src_end);

    // loop on the ghost list
    for (const auto gblock : (*ghost_list)) {
        GridBlock* ngh_block = gblock->data_src();
        // get info that change with the GP: where to put and from where to take it
        MemLayout* block_trg = gblock;
        real_p     data_src  = ngh_block->data(fid, ida_);
        // copy the information
        interp_->Copy(gblock->dlvl(), gblock->shift(), block_src, data_src, block_trg, data_trg);
    }
    //-------------------------------------------------------------------------
}

void Ghost::Compute4Block_Copy2Coarse_(const ListGBLocal* ghost_list, Field* fid, GridBlock* block_trg, real_t* ptr_trg) {
    //-------------------------------------------------------------------------
    // the source block is always the same
    const lid_t src_start[3] = {0, 0, 0};
    const lid_t src_end[3]   = {M_N, M_N, M_N};
    MemLayout*  block_src    = new SubBlock(M_GS, M_STRIDE, src_start, src_end);

    // the coarse block is computed later
    lid_t coarse_start[3];
    lid_t coarse_end[3];

    // loop on the ghost list
    for (const auto gblock : (*ghost_list)) {
        GridBlock* ngh_block = gblock->data_src();
        // set the coarse block to the correct position
        for (sid_t id = 0; id < 3; id++) {
            coarse_start[id] = CoarseFromBlock(gblock->start(id), interp_);
            coarse_end[id]   = CoarseFromBlock(gblock->end(id), interp_);
        }
        MemLayout* block_trg = new SubBlock(CoarseNGhostFront(interp_), CoarseStride(interp_), coarse_start, coarse_end);
        real_p     data_src  = ngh_block->data(fid, ida_);
        real_p     data_trg  = ptr_trg + m_zeroidx(0, block_trg);
        // interpolate, the level is 1 coarser and the shift is unchanged
        m_assert((gblock->dlvl() + 1) == 1 || (gblock->dlvl() + 1) == 0, "the difference of level MUST be 1");
        // m_log("entering interpolator with shift = %d %d %d", gblock->shift(0), gblock->shift(1), gblock->shift(2));
        // m_log("entering interpolator with srcstart = %d %d %d", block_src->start(0), block_src->start(1), block_src->start(2));
        // m_log("entering interpolator with srcend = %d %d %d", block_src->end(0), block_src->end(1), block_src->end(2));
        // m_log("entering interpolator with trgstart = %d %d %d", block_trg->start(0), block_trg->start(1), block_trg->start(2));
        // m_log("entering interpolator with trgend = %d %d %d", block_trg->end(0), block_trg->end(1), block_trg->end(2));
        interp_->Copy(gblock->dlvl() + 1, gblock->shift(), block_src, data_src, block_trg, data_trg);
    }
    //-------------------------------------------------------------------------
}
void Ghost::Compute4Block_GetRma2Myself_(const ListGBMirror* ghost_list, Field* fid, GridBlock* block_trg, real_t* data_trg) {
    //-------------------------------------------------------------------------
    // the source block is always the same
    const lid_t src_start[3] = {0, 0, 0};
    const lid_t src_end[3]   = {M_N, M_N, M_N};
    MemLayout*  block_src    = new SubBlock(M_GS, M_STRIDE, src_start, src_end);

    // loop on the ghost list
    for (const auto gblock : (*ghost_list)) {
        MPI_Aint   disp_src  = gblock->data_src();
        rank_t     disp_rank = gblock->rank();
        MemLayout* block_trg = gblock;
        // copy the information
        interp_->GetRma(gblock->dlvl(), block_src, disp_src, block_trg, data_trg, disp_rank, mirrors_window_);
    }
    //-------------------------------------------------------------------------
}

void Ghost::Compute4Block_GetRma2Coarse_(const ListGBMirror* ghost_list, Field* fid, GridBlock* block_trg, real_t* ptr_trg) {
    //-------------------------------------------------------------------------
    // the source block is always the same
    const lid_t src_start[3] = {0, 0, 0};
    const lid_t src_end[3]   = {M_N, M_N, M_N};
    MemLayout*  block_src    = new SubBlock(M_GS, M_STRIDE, src_start, src_end);

    // the coarse block is computed later
    lid_t coarse_start[3];
    lid_t coarse_end[3];

    // loop on the ghost list
    for (const auto gblock : (*ghost_list)) {
        MPI_Aint disp_src = gblock->data_src();
        rank_t disp_rank = gblock->rank();
        // set the coarse block to the correct position
        for (int id = 0; id < 3; id++) {
            coarse_start[id] = CoarseFromBlock(gblock->start(id), interp_);
            coarse_end[id]   = CoarseFromBlock(gblock->end(id), interp_);
        }
        MemLayout* block_trg = new SubBlock(CoarseNGhostFront(interp_), CoarseStride(interp_), coarse_start, coarse_end);
        real_p     data_trg  = ptr_trg + m_zeroidx(0, block_trg);

        // interpolate, the level is 1 coarser and the shift is unchanged
        m_assert((gblock->dlvl() + 1) == 1, "the difference of level MUST be 1");
        interp_->GetRma(gblock->dlvl() + 1, block_src, disp_src, block_trg, data_trg,disp_rank,mirrors_window_);
    }
    //-------------------------------------------------------------------------
}

void Ghost::Compute4Block_Myself2Coarse_(const qid_t* qid, GridBlock* cur_block, Field* fid, real_t* ptr_trg) {
    //-------------------------------------------------------------------------
    // m_log("entering myself");
    lid_t coarse_start[3];
    lid_t coarse_end[3];

    //set the coarse block to its whole domain
    for (int id = 0; id < 3; id++) {
        coarse_start[id] = 0;
        coarse_end[id]   = M_HN;
    }
    SubBlock* coarse_block = new SubBlock(CoarseNGhostFront(interp_), CoarseStride(interp_), coarse_start, coarse_end);
    // get memory details
    lid_t      shift[3]  = {0, 0, 0};
    MemLayout* block_src = cur_block;
    real_p     data_src  = cur_block->data(fid, ida_);
    MemLayout* block_trg = coarse_block;
    real_p     data_trg  = ptr_trg + m_zeroidx(0, coarse_block);
    // interpolate
    interp_->Copy(1, shift, block_src, data_src, block_trg, data_trg);

    //................................................
    // do here some physics, to completely fill the coarse block before the interpolation
    for (auto gblock : (*phys_[qid->cid])) {
        // get the direction and the corresponding bctype
        const sid_t    dir    = gblock->dir();
        const bctype_t bctype = fid->bctype(ida_, gblock->iface());
        // in the face direction, the start and the end are already correct, only the fstart changes
        lid_t fstart[3];
        coarse_start[dir] = CoarseFromBlock(gblock->start(dir), interp_);
        coarse_end[dir]   = CoarseFromBlock(gblock->end(dir), interp_);
        fstart[dir]       = CoarseFromBlock(face_start[gblock->iface()][dir], interp_);
        // in the other direction, we need to rescale the dimensions
        for (int id = 1; id < 3; id++) {
            coarse_start[(dir + id) % 3] = CoarseFromBlock(gblock->start((dir + id) % 3), interp_);
            coarse_end[(dir + id) % 3]   = CoarseFromBlock(gblock->end((dir + id) % 3), interp_);
            fstart[(dir + id) % 3]       = CoarseFromBlock(face_start[gblock->iface()][(dir + id) % 3], interp_);
        }
        // reset the coarse block and get the correct memory location
        coarse_block->Reset(CoarseNGhostFront(interp_), CoarseStride(interp_), coarse_start, coarse_end);
        real_p data_trg = ptr_trg + m_zeroidx(0, coarse_block);
        m_log("doing phyyysiiiics in dir %d from %d %d %d to %d %d %d",dir,gblock->start(0),gblock->start(1),gblock->start(2),gblock->end(0),gblock->end(1),gblock->end(2));
        m_log("doing phyyysiiiics in dir %d from %d %d %d to %d %d %d",dir,coarse_block->start(0),coarse_block->start(1),coarse_block->start(2),coarse_block->end(0),coarse_block->end(1),coarse_block->end(2));
        // get the correct face_start
        if (bctype == M_BC_EVEN) {
            EvenBoundary_4 bc = EvenBoundary_4();
            bc(gblock->iface(), fstart, cur_block->hgrid(), 0.0, coarse_block, data_trg);
        } else if (bctype == M_BC_ODD) {
            OddBoundary_4 bc = OddBoundary_4();
            bc(gblock->iface(), fstart, cur_block->hgrid(), 0.0, coarse_block, data_trg);
        } else if (bctype == M_BC_EXTRAP_3) {
            ExtrapBoundary_3 bc = ExtrapBoundary_3();
            bc(gblock->iface(), fstart, cur_block->hgrid(), 0.0, coarse_block, data_trg);
        } else if (bctype == M_BC_EXTRAP_4) {
            ExtrapBoundary_4 bc = ExtrapBoundary_4();
            bc(gblock->iface(), fstart, cur_block->hgrid(), 0.0, coarse_block, data_trg);
        } else if (bctype == M_BC_EXTRAP_5) {
            ExtrapBoundary_5 bc = ExtrapBoundary_5();
            bc(gblock->iface(), fstart, cur_block->hgrid(), 0.0, coarse_block, data_trg);
        } else if (bctype == M_BC_ZERO) {
            ZeroBoundary bc = ZeroBoundary();
            bc(gblock->iface(), fstart, cur_block->hgrid(), 0.0, coarse_block, data_trg);
        } else {
            m_assert(false, "this type of BC is not implemented yet or not valid %d", bctype);
        }
        //-------------------------------------------------------------------------
    }
}

void Ghost::Compute4Block_Refine_(const ListGBLocal* ghost_list, real_t* ptr_src, real_t* data_trg) {
    //-------------------------------------------------------------------------
    // the source block is always the same
    lid_t coarse_start[3];
    lid_t coarse_end[3];
    for (sid_t id = 0; id < 3; id++) {
        coarse_start[id] = -CoarseNGhostFront(interp_);
        coarse_end[id]   = CoarseStride(interp_) - CoarseNGhostFront(interp_);
    }
    MemLayout* block_src = new SubBlock(CoarseNGhostFront(interp_), CoarseStride(interp_), coarse_start, coarse_end);
    real_p     data_src  = ptr_src + m_zeroidx(0, block_src);
    lid_t      shift[3]  = {0, 0, 0};

    // loop on the ghost list
    for (const auto gblock : (*ghost_list)) {
        // m_log("entering interpolator with shift = %d %d %d", shift[0], shift[1], shift[2]);
        // m_log("entering interpolator with srcstart = %d %d %d", block_src->start(0), block_src->start(1), block_src->start(2));
        // m_log("entering interpolator with srcend = %d %d %d", block_src->end(0), block_src->end(1), block_src->end(2));
        // m_log("entering interpolator with trgstart = %d %d %d", gblock->start(0), gblock->start(1), gblock->start(2));
        // m_log("entering interpolator with trgend = %d %d %d", gblock->end(0), gblock->end(1), gblock->end(2));
        // interpolate, the level is 1 coarser and the shift is unchanged
        interp_->Interpolate(-1, shift, block_src, data_src, gblock, data_trg);
    }
    //-------------------------------------------------------------------------
}
void Ghost::Compute4Block_Refine_(const ListGBMirror* ghost_list, real_t* ptr_src, real_t* data_trg) {
    //-------------------------------------------------------------------------
    // the source block is always the same
    lid_t coarse_start[3];
    lid_t coarse_end[3];
    for (sid_t id = 0; id < 3; id++) {
        coarse_start[id] = -CoarseNGhostFront(interp_);
        coarse_end[id]   = CoarseStride(interp_) - CoarseNGhostFront(interp_);
    }
    MemLayout* block_src = new SubBlock(CoarseNGhostFront(interp_), CoarseStride(interp_), coarse_start, coarse_end);
    real_p     data_src  = ptr_src + m_zeroidx(0, block_src);
    lid_t      shift[3]  = {0, 0, 0};

    // loop on the ghost list
    for (const auto gblock : (*ghost_list)) {
        // interpolate, the level is 1 coarser and the shift is unchanged
        interp_->Interpolate(-1, shift, block_src, data_src, gblock, data_trg);
    }
    //-------------------------------------------------------------------------
}

// void Ghost::Compute4Block_Sibling_(const list<GhostBlock*>* ghost_list, const bool do_coarse, InterpFunction* copy, GridBlock* cur_block, Field* fid, real_t* coarse_mem) {
//     //-------------------------------------------------------------------------
//     // the source block is always the same
//     const lid_t src_start[3] = {0, 0, 0};
//     const lid_t src_end[3]   = {M_N, M_N, M_N};
//     MemLayout*  block_src    = new SubBlock(M_GS, M_STRIDE, src_start, src_end);

//     // same for the destination memory, does not change
//     real_p data_trg = cur_block->data(fid, ida_);

//     // the coarse block is computed later
//     lid_t coarse_start[3];
//     lid_t coarse_end[3];

//     // loop on the ghost list
//     for (const auto gblock : (*ghost_list)) {
//         m_assert(gblock->dlvl() == 0, "we are treating sibling, the delta level must be 0");
//         // get info that change with the GP: where to put and from where to take it
//         MemLayout* block_trg = gblock;
//         real_p     data_src  = gblock->data(fid, ida_);
//         // copy the information
//         (*copy)(0, gblock->shift(), block_src, data_src, block_trg, data_trg);

//         // we need to copy the data to the coarse representation
//         if (do_coarse) {
//             // set the coarse block to the correct position
//             for (int id = 0; id < 3; id++) {
//                 coarse_start[id] = CoarseFromBlock(gblock->start(id), interp_);
//                 coarse_end[id]   = CoarseFromBlock(gblock->end(id), interp_);
//             }
//             MemLayout* block_trg = new SubBlock(CoarseNGhostFront(interp_), CoarseStride(interp_), coarse_start, coarse_end);
//             real_p     data_trg  = coarse_mem + m_zeroidx(0, block_trg);
//             // interpolate, the level is 1 coarser and the shift is unchanged
//             m_assert((gblock->dlvl() + 1) == 1, "the difference of level MUST be 1");
//             (*copy)(gblock->dlvl() + 1, gblock->shift(), block_src, data_src, block_trg, data_trg);
//         }
//     }
//     //-------------------------------------------------------------------------
// }

// void Ghost::Compute4Block_FromParent_(const list<GhostBlock*>* ghost_list, InterpFunction* copy, GridBlock* cur_block, Field* fid, real_t* coarse_mem) {
//     //-------------------------------------------------------------------------
//     // the source block is always the same
//     const lid_t src_start[3] = {0, 0, 0};
//     const lid_t src_end[3]   = {M_N, M_N, M_N};
//     MemLayout*  block_src    = new SubBlock(M_GS, M_STRIDE, src_start, src_end);

//     // same for the destination memory, does not change
//     real_p data_trg = cur_block->data(fid, ida_);

//     // the coarse block is computed later
//     lid_t coarse_start[3];
//     lid_t coarse_end[3];

//     // loop on the ghost list
//     for (const auto gblock : (*ghost_list)) {
//         m_assert(gblock->dlvl() == 0, "we are treating sibling, the delta level must be 0");
//         // get info that change with the GP: where to put and from where to take it
//         MemLayout* block_trg = gblock;
//         real_p     data_src  = gblock->data(fid, ida_);
//         // copy the information
//         (*copy)(0, gblock->shift(), block_src, data_src, block_trg, data_trg);

//         // we need to copy the data to the coarse representation
//         if (do_coarse) {
//             // set the coarse block to the correct position
//             for (int id = 0; id < 3; id++) {
//                 coarse_start[id] = CoarseFromBlock(gblock->start(id), interp_);
//                 coarse_end[id]   = CoarseFromBlock(gblock->end(id), interp_);
//             }
//             MemLayout* block_trg = new SubBlock(CoarseNGhostFront(interp_), CoarseStride(interp_), coarse_start, coarse_end);
//             real_p     data_trg  = coarse_mem + m_zeroidx(0, block_trg);
//             // interpolate, the level is 1 coarser and the shift is unchanged
//             m_assert((gblock->dlvl() + 1) == 1, "the difference of level MUST be 1");
//             (*copy)(gblock->dlvl() + 1, gblock->shift(), block_src, data_src, block_trg, data_trg);
//         }
//     }
//     //-------------------------------------------------------------------------
// }

// void Ghost::PullFromGhost4Block_Sibling_(const qid_t *qid, GridBlock *cur_block, Field *fid,
//                                               const bool do_coarse, SubBlock *ghost_block, SubBlock *coarse_block, real_t *coarse_mem){
//     //-------------------------------------------------------------------------
//     lid_t coarse_start[3];
//     lid_t coarse_end[3];

//     // m_log("entering siblings");

//     // do the blocks on my level
//     for (auto biter = block_sibling_[qid->cid]->begin(); biter != block_sibling_[qid->cid]->end(); biter++) {
//         GhostBlock* gblock    = (*biter);
//         GridBlock*  ngh_block = gblock->block_src();
//         // memory details
//         MemLayout* block_src = ngh_block;
//         real_p     data_src  = ngh_block->data(fid, ida_);
//         MemLayout* block_trg = gblock;
//         real_p     data_trg  = cur_block->data(fid, ida_);
//         // interpolate
//         m_assert(gblock->dlvl() == 0,"we are treating sibling, the delta level must be 0");
//         interp_->Interpolate(gblock->dlvl(), gblock->shift(), block_src, data_src, block_trg, data_trg);

//         // we need to copy the data to the coarse representation
//         if (do_coarse) {
//             // m_log("doooo coarse!!");
//             // set the coarse block to the correct position
//             for (int id = 0; id < 3; id++) {
//                 coarse_start[id] = CoarseFromBlock(gblock->start(id),interp_);
//                 coarse_end[id]   = CoarseFromBlock(gblock->end(id), interp_);
//             }
//             coarse_block->Reset(CoarseNGhostFront(interp_), CoarseStride(interp_), coarse_start, coarse_end);
//             // memory details
//             MemLayout* block_src = ngh_block;
//             real_p     data_src  = ngh_block->data(fid, ida_);
//             MemLayout* block_trg = coarse_block;
//             real_p     data_trg  = coarse_mem + m_zeroidx(0, coarse_block);
//             // interpolate, the level is 1 coarser and the shift is unchanged
//             m_assert((gblock->dlvl() + 1) == 1, "the difference of level MUST be 1");
//             interp_->Copy(gblock->dlvl() + 1, gblock->shift(), block_src, data_src, block_trg, data_trg);
//         }
//     }

//     // do the ghosts, on my level
//     for (auto biter = ghost_sibling_[qid->cid]->begin(); biter != ghost_sibling_[qid->cid]->end(); biter++) {
//         GhostBlock* gblock = (*biter);
//         // memory details
//         MemLayout* block_src = ghost_block;
//         real_p     data_src  = gblock->data_src();
//         MemLayout* block_trg = gblock;
//         real_p     data_trg  = cur_block->data(fid, ida_);
//         // interpolate
//         m_assert(gblock->dlvl() == 0,"we are treating sibling, the delta level must be 0");
//         interp_->Interpolate(gblock->dlvl(), gblock->shift(), block_src, data_src, block_trg, data_trg);

//         // we need to interpolate on the coarse version of myself as well
//         if (do_coarse) {
//             // set the coarse block to the correct position
//             for (int id = 0; id < 3; id++) {
//                 coarse_start[id] = CoarseFromBlock(gblock->start(id),interp_);
//                 coarse_end[id]   = CoarseFromBlock(gblock->end(id), interp_);
//             }
//             coarse_block->Reset(CoarseNGhostFront(interp_), CoarseStride(interp_), coarse_start, coarse_end);
//             // memory details
//             MemLayout* block_src = ghost_block;
//             real_p     data_src  = gblock->data_src();
//             MemLayout* block_trg = coarse_block;
//             real_p     data_trg  = coarse_mem + m_zeroidx(0, coarse_block);
//             // interpolate, the level is 1 coarser and the shift is unchanged
//             m_assert(gblock->dlvl()==1,"the difference of level MUST be 1");
//             interp_->Interpolate(gblock->dlvl() + 1, gblock->shift(), block_src, data_src, block_trg, data_trg);
//         }
//     }
//     //-------------------------------------------------------------------------
// // }

// void Ghost::PullFromGhost4Block_FromParent_(const qid_t* qid, GridBlock* cur_block, Field* fid,
//                                            const bool do_coarse, SubBlock* ghost_block, SubBlock* coarse_block, real_t* coarse_mem) {
//     //-------------------------------------------------------------------------
//     // m_log("entering from parent");
//     lid_t coarse_start[3];
//     lid_t coarse_end[3];
//     // copy the coarse blocks to the coarse representation
//     for (auto biter = block_parent_[qid->cid]->begin(); biter != block_parent_[qid->cid]->end(); biter++) {
//         GhostBlock* gblock    = (*biter);
//         GridBlock*  ngh_block = gblock->block_src();
//         // setup the coarse sublock to the position
//         for (int id = 0; id < 3; id++) {
//             coarse_start[id] = CoarseFromBlock(gblock->start(id),interp_);
//             coarse_end[id]   = CoarseFromBlock(gblock->end(id), interp_);
//         }
//         coarse_block->Reset(CoarseNGhostFront(interp_), CoarseStride(interp_), coarse_start, coarse_end);
//         // memory details
//         MemLayout* block_src = ngh_block;
//         real_p     data_src  = ngh_block->data(fid, ida_);
//         MemLayout* block_trg = coarse_block;
//         real_p     data_trg  = coarse_mem + m_zeroidx(0, coarse_block);
//         // interpolate, the level is 1 coarser and the shift is unchanged
//         m_assert(gblock->dlvl() + 1 == 0, "the gap in level has to be 0");
//         interp_->Copy(gblock->dlvl() + 1, gblock->shift(), block_src, data_src, block_trg, data_trg);
//     }
//     // copy the ghost into the coarse representation
//     for (auto biter = ghost_parent_[qid->cid]->begin(); biter != ghost_parent_[qid->cid]->end(); biter++) {
//         m_assert(false,"cannot be here!!");
//         GhostBlock* gblock = (*biter);
//         // update the coarse subblock
//         for (int id = 0; id < 3; id++) {
//             coarse_start[id] = CoarseFromBlock(gblock->start(id), interp_);
//             coarse_end[id]   = CoarseFromBlock(gblock->end(id), interp_);
//         }
//         coarse_block->Reset(CoarseNGhostFront(interp_), CoarseStride(interp_), coarse_start, coarse_end);
//         // memory details
//         MemLayout* block_src = ghost_block;
//         real_p     data_src  = gblock->data_src();
//         MemLayout* block_trg = coarse_block;
//         real_p     data_trg  = coarse_mem + m_zeroidx(0, coarse_block);
//         // interpolate, the level is 1 coarser and the shift is unchanged
//         m_assert(gblock->dlvl() + 1 == 0, "the gap in level has to be 0");
//         interp_->Copy(gblock->dlvl() + 1, gblock->shift(), block_src, data_src, block_trg, data_trg);
//     }
//     // m_log("copy of the neighbor is done, now compute the fine info");
//     //---------------------
//     // reset the coarse block to the correct position
//     for (int id = 0; id < 3; id++) {
//         coarse_start[id] = -CoarseNGhostFront(interp_);
//         coarse_end[id]   = CoarseStride(interp_) - CoarseNGhostFront(interp_);
//     }
//     coarse_block->Reset(CoarseNGhostFront(interp_), CoarseStride(interp_), coarse_start, coarse_end);
//     // now that we have everything, we compute the fine ghosts
//     for (auto biter = block_parent_[qid->cid]->begin(); biter != block_parent_[qid->cid]->end(); biter++) {
//         GhostBlock* gblock    = (*biter);
//         lid_t       shift[3]  = {0, 0, 0};
//         MemLayout*  block_src = coarse_block;
//         real_p      data_src  = coarse_mem + m_zeroidx(0, coarse_block);
//         MemLayout*  block_trg = gblock;
//         real_p      data_trg  = cur_block->data(fid, ida_);
//         // interpolate
//         interp_->Interpolate(-1, shift, block_src, data_src, block_trg, data_trg);
//     }
//     // copy the ghost into the coarse representation
//     for (auto biter = ghost_parent_[qid->cid]->begin(); biter != ghost_parent_[qid->cid]->end(); biter++) {
//         GhostBlock* gblock    = (*biter);
//         lid_t       shift[3]  = {0, 0, 0};
//         MemLayout*  block_src = coarse_block;
//         real_p      data_src  = coarse_mem + m_zeroidx(0, coarse_block);
//         MemLayout*  block_trg = gblock;
//         real_p      data_trg  = cur_block->data(fid, ida_);
//         // interpolate
//         interp_->Interpolate(-1, shift, block_src, data_src, block_trg, data_trg);
//     }
//     //-------------------------------------------------------------------------
// }

// void Ghost::PullFromGhost4Block_Myself_(const qid_t* qid, GridBlock* cur_block, Field* fid, SubBlock* coarse_block, real_t* coarse_mem) {
//     //-------------------------------------------------------------------------
//     // m_log("entering myself");
//     lid_t coarse_start[3];
//     lid_t coarse_end[3];
//     //set the coarse block to its whole domain
//     for (int id = 0; id < 3; id++) {
//         coarse_start[id] = 0;
//         coarse_end[id]   = M_HN;
//     }
//     coarse_block->Reset(CoarseNGhostFront(interp_), CoarseStride(interp_), coarse_start, coarse_end);
//     // get memory details
//     lid_t      shift[3]  = {0, 0, 0};
//     MemLayout* block_src = cur_block;
//     real_p     data_src  = cur_block->data(fid, ida_);
//     MemLayout* block_trg = coarse_block;
//     real_p     data_trg  = coarse_mem + m_zeroidx(0, coarse_block);
//     // interpolate
//     interp_->Copy(1, shift, block_src, data_src, block_trg, data_trg);


//     // do here some physics, to completely fill the coarse block before the interpolation
//     for (auto piter = phys_[qid->cid]->begin(); piter != phys_[qid->cid]->end(); piter++) {
//         PhysBlock* gblock = (*piter);
//         // get the direction and the corresponding bctype
//         const sid_t    dir    = gblock->dir();
//         const bctype_t bctype = fid->bctype(ida_, gblock->iface());
//         // in the face direction, the start and the end are already correct, only the fstart changes
//         lid_t fstart[3];
//         coarse_start[dir] = CoarseFromBlock(gblock->start(dir), interp_);
//         coarse_end[dir]   = CoarseFromBlock(gblock->end(dir), interp_);
//         fstart[dir]       = CoarseFromBlock(face_start[gblock->iface()][dir], interp_);
//         // in the other direction, we need to rescale the dimensions
//         for (int id = 1; id < 3; id++) {
//             coarse_start[(dir + id) % 3] = CoarseFromBlock(gblock->start((dir + id) % 3), interp_);
//             coarse_end[(dir + id) % 3]   = CoarseFromBlock(gblock->end((dir + id) % 3), interp_);
//             fstart[(dir + id) % 3]       = CoarseFromBlock(face_start[gblock->iface()][(dir + id) % 3], interp_);
//         }
//         // reset the coarse block and get the correct memory location
//         coarse_block->Reset(CoarseNGhostFront(interp_), CoarseStride(interp_), coarse_start, coarse_end);
//         real_p data_trg = coarse_mem + m_zeroidx(0, coarse_block);
//         // get the correct face_start
//         if (bctype == M_BC_EVEN) {
//             EvenBoundary_4 bc = EvenBoundary_4();
//             bc(gblock->iface(), fstart, cur_block->hgrid(), 0.0, coarse_block, data_trg);
//         } else if (bctype == M_BC_ODD) {
//             OddBoundary_4 bc = OddBoundary_4();
//             bc(gblock->iface(), fstart, cur_block->hgrid(), 0.0, coarse_block, data_trg);
//         } else if (bctype == M_BC_EXTRAP_3) {
//             ExtrapBoundary_3 bc = ExtrapBoundary_3();
//             bc(gblock->iface(), fstart, cur_block->hgrid(), 0.0, coarse_block, data_trg);
//         } else if (bctype == M_BC_EXTRAP_4) {
//             ExtrapBoundary_4 bc = ExtrapBoundary_4();
//             bc(gblock->iface(), fstart, cur_block->hgrid(), 0.0, coarse_block, data_trg);
//         } else if (bctype == M_BC_EXTRAP_5) {
//             ExtrapBoundary_5 bc = ExtrapBoundary_5();
//             bc(gblock->iface(), fstart, cur_block->hgrid(), 0.0, coarse_block, data_trg);
//         } else if (bctype == M_BC_ZERO) {
//             ZeroBoundary bc = ZeroBoundary();
//             bc(gblock->iface(), fstart, cur_block->hgrid(), 0.0, coarse_block, data_trg);
//         } else {
//             m_assert(false, "this type of BC is not implemented yet or not valid %d", bctype);
//         }
//         //-------------------------------------------------------------------------
//     }
// }

// void Ghost::PullFromGhost4Block_ToParent_(const qid_t* qid, GridBlock* cur_block, Field* fid,
//                                          SubBlock* ghost_block, SubBlock* coarse_block, real_t* coarse_mem) {
//     //-------------------------------------------------------------------------
//     lid_t coarse_start[3];
//     lid_t coarse_end[3];

//     // m_log("entering to parent");

//     // copy the coarse blocks to the coarse representation
//     for (auto biter = block_children_[qid->cid]->begin(); biter != block_children_[qid->cid]->end(); biter++) {
//         GhostBlock* gblock    = (*biter);
//         GridBlock*  ngh_block = gblock->block_src();
//         // setup the coarse sublock to the position
//         // for (int id = 0; id < 3; id++) {
//         //     coarse_start[id] = CoarseFromBlock(gblock->start(id), cgs(interp_));
//         //     coarse_end[id]   = CoarseFromBlock(gblock->end(id), cgs(interp_));
//         // }
//         // coarse_block->Reset(cgs(interp_), CoarseStride(interp_), coarse_start, coarse_end);
//         // memory details
//         MemLayout* block_src = cur_block;
//         real_p     data_src  = cur_block->data(fid, ida_);
//         MemLayout* block_trg = gblock;
//         real_p     data_trg  = ngh_block->data(fid, ida_);
//         // interpolate, the level is 1 coarser and the shift is unchanged
//         interp_->Interpolate(gblock->dlvl(), gblock->shift(), block_src, data_src, block_trg, data_trg);
//     }
//     // copy the ghost into the coarse representation
//     for (auto biter = ghost_children_[qid->cid]->begin(); biter != ghost_children_[qid->cid]->end(); biter++) {
//         m_assert(false, "euh this needs an MPI layer!!");
//         // GhostBlock* gblock = (*biter);
//         // // update the coarse subblock
//         // for (int id = 0; id < 3; id++) {
//         //     coarse_start[id] = CoarseFromBlock(gblock->start(id), cgs(interp_));
//         //     coarse_end[id]   = CoarseFromBlock(gblock->end(id), cgs(interp_));
//         // }
//         // coarse_block->Reset(cgs(interp_), CoarseStride(interp_), coarse_start, coarse_end);
//         // // memory details
//         // MemLayout* block_src = ghost_block;
//         // real_p     data_src  = gblock->data_src();
//         // MemLayout* block_trg = coarse_block;
//         // real_p     data_trg  = coarse_mem + m_zeroidx(0, coarse_block);
//         // // interpolate, the level is 1 coarser and the shift is unchanged
//         // m_assert(gblock->dlvl() + 1 == 0, "the gap in level has to be 0");
//         // interp_->Interpolate(gblock->dlvl() + 1, gblock->shift(), block_src, data_src, block_trg, data_trg);
//     }
//     //-------------------------------------------------------------------------
// }

// void Ghost::PullFromGhost4Block_Physics(PhysBlock* gblock, Field* fid, real_t hgrid[3], real_t* data) {
//     //-------------------------------------------------------------------------
//     // get the direction and the corresponding bctype
//     bctype_t bctype = fid->bctype(ida_, gblock->iface());
//     if (bctype == M_BC_EVEN) {
//         EvenBoundary_4 bc = EvenBoundary_4();
//         bc(gblock->iface(), hgrid, 0.0, gblock, data);
//     } else if (bctype == M_BC_ODD) {
//         OddBoundary_4 bc = OddBoundary_4();
//         bc(gblock->iface(), hgrid, 0.0, gblock, data);
//     } else if (bctype == M_BC_EXTRAP_3) {
//         ExtrapBoundary_3 bc = ExtrapBoundary_3();
//         bc(gblock->iface(), hgrid, 0.0, gblock, data);
//     } else if (bctype == M_BC_EXTRAP_4) {
//         ExtrapBoundary_4 bc = ExtrapBoundary_4();
//         bc(gblock->iface(), hgrid, 0.0, gblock, data);
//     } else if (bctype == M_BC_EXTRAP_5) {
//         ExtrapBoundary_5 bc = ExtrapBoundary_5();
//         bc(gblock->iface(), hgrid, 0.0, gblock, data);
//     } else if (bctype == M_BC_ZERO) {
//         ZeroBoundary bc = ZeroBoundary();
//         bc(gblock->iface(), hgrid, 0.0, gblock, data);
//     } else {
//         m_assert(false, "this type of BC is not implemented yet");
//     }
//     //-------------------------------------------------------------------------
// }

/**
 * @brief Loop on the blocks that are mirrors and call a gop_t operation on them
 * 
 * @param op 
 * @param field 
 */
void Ghost::LoopOnMirrorBlock_(const gop_t op, Field* field) {
    m_begin;
    m_assert(grid_->is_mesh_valid(), "mesh is not valid, unable to process");
    //-------------------------------------------------------------------------
    // get the grid info
    p8est_t*       forest = grid_->forest();
    p8est_ghost_t* ghost  = grid_->ghost();
    // const lid_t    nqlocal = ghost->mirrors.elem_count;  //number of ghost blocks

#pragma omp parallel for
    for (lid_t bid = 0; bid < n_mirror_to_send_; bid++) {
        // get the mirror quad, this is an empty quad (just a piggy3 struct)
        // p8est_quadrant_t* mirror = p8est_quadrant_array_index(&ghost->mirrors, local_to_mirrors[bid]);
        p8est_quadrant_t* mirror = p8est_quadrant_array_index(&ghost->mirrors, bid);
        p8est_tree_t*     tree   = p8est_tree_array_index(forest->trees, mirror->p.piggy3.which_tree);
        // get the id, in this case the cummulative id = mirror id
        qid_t myid;
        myid.cid = bid;
        myid.qid = mirror->p.piggy3.local_num - tree->quadrants_offset;
        myid.tid = mirror->p.piggy3.which_tree;
        // use it to retreive the actual quadrant in the correct tree
        p8est_quadrant_t* quad  = p8est_quadrant_array_index(&tree->quadrants, myid.qid);
        GridBlock*        block = *(reinterpret_cast<GridBlock**>(quad->p.user_data));
        // send the task
        (this->*op)(&myid, block, field);
    }
    //-------------------------------------------------------------------------
    m_end;
}
