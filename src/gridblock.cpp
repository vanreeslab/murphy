#include "gridblock.hpp"

#include <p8est_bits.h>

#include <algorithm>

#include "p8est_iterate.h"
#include "toolsp4est.hpp"

using std::string;
using std::unordered_map;

using GBLocal      = GhostBlock<GridBlock*>;
using GBMirror     = GhostBlock<MPI_Aint>;
using GBPhysic     = PhysBlock;
using ListGBLocal  = std::list<GBLocal*>;
using ListGBMirror = std::list<GBMirror*>;
using listGBPhysic = std::list<GBPhysic*>;

#define M_NNEIGHBORS 26

static lid_t face_start[6][3] = {{0, 0, 0}, {M_N, 0, 0}, {0, 0, 0}, {0, M_N, 0}, {0, 0, 0}, {0, 0, M_N}};

/**
 * @brief Given a face, edge or a corner returns the outgoing normal (sign)
 * 
 * @param ibidule for a face (`ibidule<6`), an edge (`6<= ibidule < 18`) or a corner (`18<= ibidule < 26`)
 * @param sign the sign of the outgoing normal
 */
static void GhostGetSign(const iface_t ibidule, real_t sign[3]) {
    // we need to find the sign = the direction of the normal:
    sign[0] = 0.0;
    sign[1] = 0.0;
    sign[2] = 0.0;

    // check depending on the plane, the edge of the corner
    if (ibidule < 6) {
        iface_t dir = ibidule / 2;
        sign[dir]   = ((ibidule % 2) == 1) ? 1.0 : -1.0;
        m_assert(fabs(sign[0]) + fabs(sign[1]) + fabs(sign[2]) == 1, "we cannot have more than 1 nonzero sign: %f %f %f (bidule=%d)", sign[0], sign[1], sign[2], ibidule);
    } else if (ibidule < 18) {
        iface_t iedge = ibidule - 6;
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
        iface_t dir  = iedge / 4;           // this is the direction of the edge
        iface_t dir1 = (dir == 0) ? 1 : 0;  // dir1 in the plane: dir1 = x if dir = y or z, or y if dir = x
        iface_t dir2 = (dir == 2) ? 1 : 2;  // dir2 in the plane: dir2 = y if dir=z, = z if dir=x or dir = y
        sign[dir1]   = ((iedge % 4) % 2) == 1 ? +1.0 : -1.0;
        sign[dir2]   = ((iedge % 4) / 2) == 1 ? +1.0 : -1.0;
        m_assert(fabs(sign[0]) + fabs(sign[1]) + fabs(sign[2]) == 2, "we cannot have more than 1 nonzero sign: %f %f %f (bidule=%d, dir = %d)", sign[0], sign[1], sign[2], ibidule, dir);
    } else {
        iface_t icorner = ibidule - 18;
        sign[0]         = (icorner % 2) == 1 ? +1.0 : -1.0;
        sign[1]         = ((icorner % 4) / 2) == 1 ? +1.0 : -1.0;
        sign[2]         = (icorner / 4) == 1 ? +1.0 : -1.0;
        m_assert(fabs(sign[0]) + fabs(sign[1]) + fabs(sign[2]) == 3, "we cannot have more than 1 nonzero sign: %f %f %f (bidule=%d)", sign[0], sign[1], sign[2], ibidule);
    }

    m_assert(sign[0] == 0.0 || sign[0] == 1.0 || sign[0] == -1.0, "wrong sign value: %e", sign[0]);
    m_assert(sign[1] == 0.0 || sign[1] == 1.0 || sign[1] == -1.0, "wrong sign value: %e", sign[1]);
    m_assert(sign[2] == 0.0 || sign[2] == 1.0 || sign[2] == -1.0, "wrong sign value: %e", sign[2]);
}

/**
 * @brief constructs a new Block given a 3D length and a position
 * 
 * @param length the length of the current block
 * @param xyz the position of the origin, i.e. the left,bottom corner, (x,y,z)
 * @param level the level of the block
 */
GridBlock::GridBlock(const real_t length, const real_t xyz[3], const sid_t level) {
    m_begin;
    //-------------------------------------------------------------------------
    level_      = level;
    status_lvl_ = 0;
    for (lda_t id = 0; id < 3; id++) {
        xyz_[id]   = xyz[id];
        hgrid_[id] = length / M_N;  // the grid spacing is still given by L/N
    }
    // init the dependencies
    n_dependency_active_ = 0;
    for (sid_t id = 0; id < P8EST_CHILDREN; ++id) {
        dependency_[id] = nullptr;
    }
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief Destroy the Grid Block
 * 
 * delete the ghost ptr if present and delete the associated fields
 * 
 */
GridBlock::~GridBlock() {
    //-------------------------------------------------------------------------
    m_verb("destruct gridBlock");
    // free the tmp memory
    if (coarse_ptr_ != nullptr) {
        m_free(coarse_ptr_);
        coarse_ptr_ = nullptr;
        m_verb("freeing the ghost ptr @ %p", coarse_ptr_);
    }

    // need to delete the fields not deleted yet
    for (auto iter = data_map_.begin(); iter != data_map_.end(); iter++) {
        m_free(iter->second);
    }
    data_map_.clear();
    //-------------------------------------------------------------------------
}

/**
 * @brief allocate the ptr for the ghost if not already exiting
 * 
 * @param memsize the memory size (in BYTES)
 */
void GridBlock::AllocateCoarsePtr(const size_t memsize) {
    //-------------------------------------------------------------------------
    if (coarse_ptr_ == nullptr) {
        coarse_ptr_ = reinterpret_cast<mem_ptr>(m_calloc(memsize));
        m_verb("allocating the ghost ptr of size %ld, @ %p", memsize, coarse_ptr_);
        m_assert(coarse_ptr_ != nullptr, "the pointer is still null, not possible");
    } else {
        m_verb("no ghost ptr allocated, already present");
    }
    //-------------------------------------------------------------------------
}

/**
 * @brief returns the (aligned!) pointer for write access that corresponds to the first point in the block, i.e. (0,0,0), for the first dimension.
 * You must use either @ref m_sidx, @ref m_midx or @ref m_idx to access any point in the memory
 * 
 * @warning this is not the same pointer as the memory pointers, because the ghost blocks are considered as negative numbers, see @ref MemLayout
 * 
 * @param fid 
 * @return real_p 
 */
data_ptr GridBlock::data(const Field* fid) {
    //-------------------------------------------------------------------------
#ifndef NDEBUG
    // check the field validity
    const auto it = data_map_.find(fid->name());
    m_assert(it != data_map_.end(), "the field %s does not exist in this block", fid->name().c_str());
    // check the alignment in memory
    data_ptr data = it->second + m_zeroidx(0, this);
    m_assert(m_isaligned(data), "M_GS and M_N have to be chosen so that (0,0,0) is aligned in memory");
    return data;
#else
    return data_map_.at(fid->name()) + m_zeroidx(0, this);
#endif
    //-------------------------------------------------------------------------
}

/**
 * @brief returns the (aligned!) pointer for write access that corresponds to the first point in the block, i.e. (0,0,0), for the given dimension.
 * You must use either @ref m_sidx, @ref m_midx or @ref m_idx to access any point in the memory
 * 
 * @warning this is not the same pointer as the memory pointers, because the ghost blocks are considered as negative numbers, see @ref MemLayout
 * 
 * @param fid the field
 * @param ida the required dimension
 * @return real_p the memory adress, we ensure its alignement
 */
data_ptr GridBlock::data(const Field* fid, const sid_t ida) {
    //-------------------------------------------------------------------------
#ifndef NDEBUG
    // check the field validity
    auto it = data_map_.find(fid->name());
    m_assert(it != data_map_.end(), "the field \"%s\" does not exist in this block", fid->name().c_str());
    // check the alignment in memory
    data_ptr data = it->second + m_zeroidx(ida, this);
    m_assert(m_isaligned(data), "M_GS = %d and M_N = %d have to be chosen so that (0,0,0) is aligned in memory: ida = %d -> o", M_GS, M_N, ida);
    return data;
#else
    return data_map_.at(fid->name()) + m_zeroidx(ida, this);
#endif
    //-------------------------------------------------------------------------
}

/**
 * @brief returns the (aligned!) pointer for write access that corresponds to the actual data pointer, for the first dimension.
 * 
 * @warning do not confuse with @ref data() functions
 * 
 * @param fid 
 * @return real_p 
 */
mem_ptr GridBlock::pointer(const Field* fid) {
#ifndef NDEBUG
    // check the field validity
    auto it = data_map_.find(fid->name());
    m_assert(it != data_map_.end(), "the field %s does not exist in this block", fid->name().c_str());
    // check the alignment in memory
    mem_ptr data = it->second;
    m_assert(m_isaligned(data), "M_GS and M_N have to be chosen so that (0,0,0) is aligned in memory");
    return data;
#else
    return data_map_.at(fid->name());
#endif
}

/**
 * @brief returns the (aligned!) pointer for write access that corresponds to the actual data pointer, for the first dimension.
 * 
 * @warning do not confuse with @ref data() functions
 * 
 * @param fid the field
 * @param ida the required dimension
 * @return real_p the memory adress
 */
mem_ptr GridBlock::pointer(const Field* fid, const sid_t ida) {
#ifndef NDEBUG
    // check the field validity
    auto it = data_map_.find(fid->name());
    m_assert(it != data_map_.end(), "the field \"%s\" does not exist in this block", fid->name().c_str());
    // check the alignment in memory
    mem_ptr data = it->second + m_blockmemsize(ida);
    m_assert(m_isaligned(data), "M_GS = %d and M_N = %d have to be chosen so that (0,0,0) is aligned in memory: ida = %d -> o", M_GS, M_N, ida);
    return data;
#else
    return data_map_.at(fid->name()) + m_blockmemsize(ida);
#endif
}

/**
 * @brief adds a field to the block if it doesn't exist already
 * 
 * @param fid 
 */
void GridBlock::AddField(Field* fid) {
    //-------------------------------------------------------------------------
    string name = fid->name();
    // try to find the field
    auto it = data_map_.find(name);
    // if not found, create it
    if (it == data_map_.end()) {
        m_verb("adding field %s to the block (dim=%d)", name.c_str(),fid->lda());
        data_map_[name] = (real_p)m_calloc(m_blockmemsize(fid->lda()) * sizeof(real_t));
    } else {
        m_verb("field %s already in the block (dim=%d)", name.c_str(),fid->lda());
    }
    //-------------------------------------------------------------------------
}

/**
 * @brief add all the fields contained in the map to the current block, if they do not exist already
 * 
 * @param fields 
 */
void GridBlock::AddFields(const unordered_map<string, Field*>* fields) {
    //-------------------------------------------------------------------------
    // remember if I need to free the memory:
    for (auto iter = fields->begin(); iter != fields->end(); iter++) {
        AddField(iter->second);
    }
    //-------------------------------------------------------------------------
}

/**
 * @brief remove the field from the current block if it exists
 * 
 * @param fid the field to remove
 */
void GridBlock::DeleteField(Field* fid) {
    //-------------------------------------------------------------------------
    string name = fid->name();
    // try to find the field
    auto it = data_map_.find(name);
    // if not found, create it
    if (it != data_map_.end()) {
        m_verb("deleting field %s to the block", name.c_str());
        m_free(data_map_[name]);
        data_map_.erase(name);
    } else {
        m_verb("no field %s in the block", name.c_str());
    }
    //-------------------------------------------------------------------------
}

/**
 * @brief compute the criterion and update GridBlock::status_lvl_ if we need to refine/coarsen
 * 
 * The criterion of the block must always be:
 *  rtol >= criterion >= ctol
 * 
 * if citerion > rtol, we must refine, whatever the other field's dimension value
 * if citerion < ctol, we can coarsen if this is satisfy by every dimension of the field
 * 
 * We loop over the field dimensions and we take the conservative approach:
 * if we need to refine, we always do it, even if other dimension may be coarsened
 * if we need to coarsen, we make sure that every dimension feels ok with it.
 * 
 * @param interp the wavelet to use to compute the criterion
 * @param rtol the refinement tolerance
 * @param ctol the coarsening tolerance
 * @param field_citerion the field that will be used as a criterion
 */
void GridBlock::UpdateStatusCriterion(const Wavelet* interp, const real_t rtol, const real_t ctol, const Field* field_citerion, Prof* profiler) {
    m_assert(rtol > ctol, "the refinement tolerance must be > the coarsening tolerance: %e vs %e", rtol, ctol);
    m_assert(status_lvl_ == 0, "trying to update a status which is already updated");
    m_assert(field_citerion->ghost_status(),"the ghost of <%s> must be up-to-date",field_citerion->name().c_str());
    //-------------------------------------------------------------------------
    m_profStart(profiler, "adapt detail");

    // I need to visit every dimension and determine if we have to refine and/or coarsen.
    // afterthat we choose given the conservative approach
    bool coarsen = true;
    for (lda_t ida = 0; ida < field_citerion->lda(); ida++) {
        // go to the computation
        data_ptr data = this->data(field_citerion, ida);
        real_t   norm = interp->Criterion(this, data);

        // if the norm is bigger than the refinement tol, we must refine
        bool refine = norm > rtol;
        if (refine) {
            status_lvl_ = +1;
            // finito
            m_profStop(profiler, "adapt detail");
            return;
        }
        // if one dimension is preventing the coarsening, register
        coarsen &= (norm < ctol);
    }
    // if every field is ok to be coarsened, i.e. the coarsen bool is still true after everything, we coarsen
    status_lvl_ = (coarsen) ? -1 : 0;

    m_profStop(profiler, "adapt detail");
    //-------------------------------------------------------------------------
}

/**
 * @brief Compute the detail coefficient and store them in the field details
 * 
 * @param interp the wavelet object
 * @param criterion the criterion field
 * @param details the detail field with the compute detail values
 */
void GridBlock::ComputeDetails(const Wavelet* interp, const Field* criterion, const Field* details) {
    m_assert(criterion->lda() == details->lda(), "field <%s> and <%s> must have the same size", criterion->name().c_str(), details->name().c_str());
    //-------------------------------------------------------------------------
    for (lda_t ida = 0; ida < criterion->lda(); ida++) {
        interp->WriteDetails(this, this->data(criterion, ida), this->data(details, ida));
    }
    //-------------------------------------------------------------------------
}

/**
 * @brief resolve the dependency list created while adapting the mesh (see cback_UpdateDependency() ) by interpolating the needed blocks
 * 
 * @param interp the Wavelet object to use for the interpolation/refinement
 * @param field_start the first field to take into account
 * @param field_end the last field to take into account
 */
void GridBlock::SolveDependency(const Wavelet* interp, std::unordered_map<std::string, Field*>::const_iterator field_start, std::unordered_map<std::string, Field*>::const_iterator field_end, Prof* profiler) {
    m_assert(n_dependency_active_ == 0 || n_dependency_active_ == 1 || n_dependency_active_ == P8EST_CHILDREN, "wrong value for n_dependency_active_");
    m_profStart(profiler, "adapt dependency");
    //-------------------------------------------------------------------------
    if (n_dependency_active_ == 1) {
        // if I get only one dependency, I am a child and I need refinement from my parent
        // we are the children of a group, we go to the parent
        GridBlock* root = this->PopDependency(0);
        m_assert(n_dependency_active_ == 0, "I should be empty now");

        // from the parent, we interpolate to me
        int childid = p4est_GetChildID(xyz_, level_);

        // get the shift given the child id
        const lid_t shift[3]     = {M_HN * ((childid % 2)), M_HN * ((childid % 4) / 2), M_HN * ((childid / 4))};
        const lid_t src_start[3] = {shift[0] - M_GS, shift[1] - M_GS, shift[2] - M_GS};
        const lid_t src_end[3]   = {shift[0] + M_HN + M_GS, shift[1] + M_HN + M_GS, shift[2] + M_HN + M_GS};
        SubBlock    mem_src(M_GS, M_STRIDE, src_start, src_end);

        // for every field on my parent, interpolate it
        m_assert(data_map_.size() == 0, "the block should be empty here");
        for (auto fid = field_start; fid != field_end; ++fid) {
            auto current_field = fid->second;
            // allocate the field on me (has not been done yet)
            this->AddField(current_field);
            // refine for every dimension
            if (!current_field->is_temp()) {
                for (sid_t ida = 0; ida < current_field->lda(); ida++) {
                    // get the pointers
                    interp->Interpolate(-1, shift, &mem_src, root->data(current_field, ida), this, this->data(current_field, ida));
                }
            }
        }
        // remove my ref from the parent
        GridBlock* this_should_be_me = root->PopDependency(childid);
        m_assert(this_should_be_me == this, "this should be me");

        // destroy my parent if I was the last one
        if (root->n_dependency_active() == 0) {
            delete (root);
        }
    } else if (n_dependency_active_ == P8EST_CHILDREN) {
        // I have 8 deps, I am a root, waiting data from coarsening of my children

        //allocate the new fields
        m_assert(data_map_.size() == 0, "the block should be empty here");
        for (auto fid = field_start; fid != field_end; ++fid) {
            auto current_field = fid->second;
            // allocate the field on me (has not been done yet)
            this->AddField(current_field);
        }

        // I am a parent and I need to fillout my children
        for (sid_t childid = 0; childid < P8EST_CHILDREN; ++childid) {
            GridBlock* child_block = this->PopDependency(childid);
            m_assert(child_block->level() - this->level() == 1, "the child block is not a child");
            m_assert(childid == p4est_GetChildID(child_block->xyz(), child_block->level()), "the two ids must match");

            // get the shift for me
            const lid_t shift[3]     = {-M_N * ((childid % 2)), -M_N * ((childid % 4) / 2), -M_N * ((childid / 4))};
            const lid_t trg_start[3] = {M_HN * ((childid % 2)), M_HN * ((childid % 4) / 2), M_HN * ((childid / 4))};
            const lid_t trg_end[3]   = {trg_start[0] + M_HN, trg_start[1] + M_HN, trg_start[2] + M_HN};
            SubBlock    mem_trg(M_GS, M_STRIDE, trg_start, trg_end);
            // and an extended source block for my child
            const lid_t src_start[3] = {-M_GS, -M_GS, -M_GS};
            const lid_t src_end[3]   = {M_N + M_GS, M_N + M_GS, M_N + M_GS};
            SubBlock    mem_src(M_GS, M_STRIDE, src_start, src_end);

            // for every field, we interpolate it
            for (auto fid = field_start; fid != field_end; ++fid) {
                auto current_field = fid->second;
                if (!current_field->is_temp()) {
                    // interpolate for every dimension
                    for (sid_t ida = 0; ida < current_field->lda(); ida++) {
                        // get the pointers
                        interp->Interpolate(1, shift, &mem_src, child_block->data(current_field, ida), &mem_trg, this->data(current_field, ida));
                    }
                }
            }
            // remove the ref to the child and my ref in the child
            child_block->PopDependency(0);
            m_assert(child_block->n_dependency_active() == 0, "the child must be empty now");
            // delete the block as we are the only one to access it
            delete (child_block);
        }
    }
    m_profStop(profiler, "adapt dependency");
    //-------------------------------------------------------------------------
}

/**
 * @brief remove the dependency for the given child
 * 
 * @param child_id the id of the dependency to remove, 0 if the parent is removed
 * @return GridBlock* 
 */
GridBlock* GridBlock::PopDependency(const sid_t child_id) {
    m_assert(0 <= child_id && child_id < P8EST_CHILDREN, "child id is out of bound");
    m_assert(dependency_[child_id] != nullptr, "there is nobody here: dep[%d] = %p", child_id, dependency_[child_id]);
    //-------------------------------------------------------------------------
    --n_dependency_active_;
    GridBlock* block      = dependency_[child_id];
    dependency_[child_id] = nullptr;
    return block;
    //-------------------------------------------------------------------------
}

/**
 * @brief add a dependency for a given child
 * 
 * @param child_id the id of the dependency to add, 0 if the dependent is the parent
 * @param dependent_block the pointer to the dependent block
 */
void GridBlock::PushDependency(const sid_t child_id, GridBlock* dependent_block) {
    m_assert(0 <= child_id && child_id < P8EST_CHILDREN, "child id is out of bound");
    m_assert(dependency_[child_id] == nullptr, "there is already someone here");
    //-------------------------------------------------------------------------
    ++n_dependency_active_;
    dependency_[child_id] = dependent_block;
    //-------------------------------------------------------------------------
}

/**
 * @brief build the list of ghosts
 * 
 * @param qid the current quadrant ID
 * @param grid the grid to use to recover the ghosts etc
 * @param interp the wavelet to use (for coarse ghost sizes etc)
 * @param local2disp_window the displacement information for RMA
 */
void GridBlock::GhostInitLists(const qid_t* qid, const ForestGrid* grid, const Wavelet* interp, MPI_Win local2disp_window) {
    //-------------------------------------------------------------------------
    // allocate the ghost pointer
    AllocateCoarsePtr(interp->CoarseMemSize());

    //................................................
    p8est_t*              forest  = grid->p4est_forest();
    p8est_mesh_t*         mesh    = grid->p4est_mesh();
    p8est_ghost_t*        ghost   = grid->p4est_ghost();
    p8est_connectivity_t* connect = forest->connectivity;

    std::list<qdrt_t*> ngh_list;
    std::list<rank_t>  rank_list;

    //................................................
    // get the number of ghost and the min/max of a block
    lid_t  block_min[3], block_max[3];
    real_t block_len[3];
    real_t coarse_hgrid[3];
    for (lda_t id = 0; id < 3; id++) {
        // set the number of ghost to compute
        block_min[id]    = -interp->nghost_front();
        block_max[id]    = M_N + interp->nghost_back();
        block_len[id]    = p4est_QuadLen(level());
        coarse_hgrid[id] = p4est_QuadLen(level()) / M_HN;
    }

    rank_t my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    for (iface_t ibidule = 0; ibidule < M_NNEIGHBORS; ibidule++) {
        //................................................
        p4est_GetNeighbor(forest, connect, ghost, mesh, qid->tid, qid->qid, ibidule, &ngh_list, &rank_list);
        const iblock_t nghosts = ngh_list.size();

        //................................................
        // no ghosts? then is a physical BC
        if (nghosts == 0) {
            // we only apply the physics to entire faces
            if (ibidule < 6) {
                PhysBlock* pb = new PhysBlock(ibidule, this, interp->nghost_front(), interp->nghost_back());
                //#pragma omp critical
                phys_.push_back(pb);
                m_verb("I found a physical boundary ghost!\n");
            }
            // else, the edges and corners will be filled through the face
        }

        //................................................
        // this is a real block or a ghost
        for (iblock_t nid = 0; nid < nghosts; ++nid) {
            qdrt_t*    nghq     = ngh_list.back();
            rank_t     ngh_rank = rank_list.back();
            const bool isghost  = (ngh_rank != my_rank);
            m_verb("reading th list: adress: %p  and rank %d -> is ghost? %d", nghq, ngh_rank, isghost);

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
                p8est_qcoord_to_vertex(grid->p4est_connect(), nghq->p.piggy3.which_tree, nghq->x, nghq->y, nghq->z, ngh_pos);
            }
            // fix the shift in coordinates needed if the domain is periodic
            for (lda_t id = 0; id < 3; ++id) {
                // if we are periodic, we overwrite the position in the direction of the normal !!ONLY!!
                // since it is my neighbor in this normal direction, I am 100% sure that it's origin corresponds to the end of my block
                const real_t to_replace = sign[id] * sign[id] * grid->domain_periodic(id);  // is (+-1)^2 = +1 if we need to replace it, 0.0 otherwize
                // get the expected position
                const real_t expected_pos = xyz(id) + (sign[id] > 0.5) * p4est_QuadLen(level()) - (sign[id] < -0.5) * p4est_QuadLen(nghq->level);
                // we override the position if a replacement is needed only
                ngh_pos[id] = to_replace * expected_pos + (1.0 - to_replace) * ngh_pos[id];
            }
            // get the hgrid
            const real_t ngh_len[3]   = {p4est_QuadLen(nghq->level), p4est_QuadLen(nghq->level), p4est_QuadLen(nghq->level)};
            const real_t ngh_hgrid[3] = {p4est_QuadLen(nghq->level) / M_N, p4est_QuadLen(nghq->level) / M_N, p4est_QuadLen(nghq->level) / M_N};

            //................................................
            // create the new block and push back
            if (!isghost) {
                // associate the corresponding neighboring block
                GridBlock* ngh_block = *(reinterpret_cast<GridBlock**>(nghq->p.user_data));

                // register the gb in a list
                if (nghq->level == level()) {
                    // sibling: source = neighbor GridBlock, target = me
                    GBLocal* gb = new GBLocal(/* source */ ngh_block->level(), ngh_pos, ngh_hgrid, ngh_len,
                                              /* traget */ level_, xyz_, hgrid_, block_min, block_max, gs(), stride(), -1);
                    gb->data_src(ngh_block);
                    //#pragma omp critical
                    local_sibling_.push_back(gb);
                } else if (nghq->level < level()) {
                    // parent: source = neighbor, target = me
                    GBLocal* gb = new GBLocal(/* source */ ngh_block->level(), ngh_pos, ngh_hgrid, ngh_len,
                                              /* target */ level_, xyz_, hgrid_, block_min, block_max, gs(), stride(), -1);
                    gb->data_src(ngh_block);
                    //#pragma omp critical
                    local_parent_.push_back(gb);
                    // the children: the source = the coarse myself, target = my neighbor
                    GBLocal* invert_gb = new GBLocal(/* source */ level_ - 1, xyz_, coarse_hgrid, block_len,
                                                     /* target */ ngh_block->level(), ngh_pos, ngh_hgrid, block_min, block_max, gs(), stride(), -1);
                    invert_gb->data_src(ngh_block);
                    //#pragma omp critical
                    local_parent_reverse_.push_back(invert_gb);

                } else {
                    m_assert((nghq->level - level_) == 1, "The delta level is not correct: %d - %d", nghq->level, level_);
                }
            }
            //................................................
            else {
                // get the local number in the remote rank and the remote rank
                rank_t ngh_local_id = nghq->p.piggy3.local_num;
                // rank_t ngh_rank     = p4est_GetOwnerFromGhost(forest, nghq);
                m_assert(ngh_rank > -1, "p4est unable to recover the rank... baaaad news");

                // register the ghost block in a list
                //................................................
                if (nghq->level == level_) {
                    // create the new mirror block
                    // GBMirror* gb = new GBMirror(block->level(), block->xyz(), block->hgrid(), nghq->level, ngh_pos, ngh_block->hgrid(), block_len, block_min, block_max, block->gs(), block->stride(), ngh_rank);
                    // sibling: source = neighbor GridBlock, target = me
                    GBMirror* gb = new GBMirror(/* source */ nghq->level, ngh_pos, ngh_hgrid, ngh_len,
                                                /* target */ level(), xyz(), hgrid(), block_min, block_max, gs(), stride(), ngh_rank);
                    // ask the displacement (will be available later, when completing the call
                    MPI_Get(gb->data_src_ptr(), 1, MPI_AINT, ngh_rank, ngh_local_id, 1, MPI_AINT, local2disp_window);
                    //#pragma omp critical
                    ghost_sibling_.push_back(gb);
                }
                //................................................
                else if (nghq->level < level_) {
                    // TODO: change this to the coarse block to avoid any cast afterwards...
                    // parent: source = neighbor, target = me
                    GBMirror* gb = new GBMirror(/* source */ nghq->level, ngh_pos, ngh_hgrid, ngh_len,
                                                /* target */ level(), xyz(), hgrid(), block_min, block_max, gs(), stride(), ngh_rank);
                    // ask the displacement (will be available later, when completing the call
                    MPI_Get(gb->data_src_ptr(), 1, MPI_AINT, ngh_rank, ngh_local_id, 1, MPI_AINT, local2disp_window);
                    //#pragma omp critical
                    ghost_parent_.push_back(gb);
                    // I compute my own contribution to my neighbor ghost points
                    GBMirror* invert_gb = new GBMirror(/* source */ level() - 1, xyz(), coarse_hgrid, block_len,
                                                       /* target */ nghq->level, ngh_pos, ngh_hgrid, block_min, block_max, gs(), stride(), ngh_rank);
                    MPI_Get(invert_gb->data_src_ptr(), 1, MPI_AINT, ngh_rank, ngh_local_id, 1, MPI_AINT, local2disp_window);
                    //#pragma omp critical
                    ghost_parent_reverse_.push_back(invert_gb);
                }
                //................................................
                else if (nghq->level > level_) {
                    const real_t ngh_hgrid_coarse[3] = {ngh_len[0] / M_HN, ngh_len[1] / M_HN, ngh_len[2] / M_HN};
                    // children: source = coarse version of my neighbor, target = myself
                    GBMirror* gb = new GBMirror(/* source */ nghq->level - 1, ngh_pos, ngh_hgrid_coarse, ngh_len,
                                                /* target */ level(), xyz(), hgrid(), block_min, block_max, gs(), stride(), ngh_rank);
                    //#pragma omp critical
                    ghost_children_.push_back(gb);
                }
                //................................................
                else {
                    m_assert(false, "The delta level is not correct: %d - %d", nghq->level, level());
                }
            }
            // pop the last element, will go to the next one
            ngh_list.pop_back();
            rank_list.pop_back();
        }
    }
    //-------------------------------------------------------------------------
}

/**
 * @brief free and clear the ghost list
 * 
 */
void GridBlock::GhostFreeLists() {
    //-------------------------------------------------------------------------
    // clear the ghost lists
    auto remove_block = [](auto block) { delete (block); };
    std::for_each(local_sibling_.begin(), local_sibling_.end(), remove_block);
    std::for_each(local_parent_.begin(), local_parent_.end(), remove_block);
    std::for_each(local_parent_reverse_.begin(), local_parent_reverse_.end(), remove_block);
    std::for_each(ghost_sibling_.begin(), ghost_sibling_.end(), remove_block);
    std::for_each(ghost_parent_.begin(), ghost_parent_.end(), remove_block);
    std::for_each(ghost_children_.begin(), ghost_children_.end(), remove_block);
    std::for_each(ghost_parent_reverse_.begin(), ghost_parent_reverse_.end(), remove_block);
    std::for_each(phys_.begin(), phys_.end(), remove_block);

    // clear the lists
    local_sibling_.clear();
    local_parent_.clear();
    local_parent_reverse_.clear();
    ghost_sibling_.clear();
    ghost_parent_.clear();
    ghost_children_.clear();
    ghost_parent_reverse_.clear();
    phys_.clear();
    //-------------------------------------------------------------------------
}

/**
 * @brief Do the first part of GhostGet and start the RMA requests
 * 
 * - get the siblings values (both local and RMA)
 * - get the coarser values to my temp array
 * 
 * @param field the field to interpolate
 * @param ida the current dimension
 * @param interp the wavelet
 * @param mirrors_window the window where to find the mirrors
 */
void GridBlock::GhostGet_Post(const Field* field, const lda_t ida, const Wavelet* interp, MPI_Win mirrors_window) {
    //-------------------------------------------------------------------------
    // get the sibligngs
    {
        const SubBlock bsrc_neighbor(M_GS, M_STRIDE, 0, M_N);
        const data_ptr data_trg = data(field, ida);

        // start the RMA for the sibling neighbors
        for (const auto gblock : ghost_sibling_) {
            const MPI_Aint   disp_src  = gblock->data_src();
            const rank_t     disp_rank = gblock->rank();
            const MemLayout* block_trg = gblock;
            // copy the information
            interp->GetRma(gblock->dlvl(), gblock->shift(), &bsrc_neighbor, disp_src, block_trg, data_trg, disp_rank, mirrors_window);
        }

        // get the copy of local siblings
        for (const auto gblock : (local_sibling_)) {
            GridBlock*       ngh_block = gblock->data_src();
            const data_ptr   data_src  = ngh_block->data(field, ida);
            const MemLayout* block_trg = gblock;
            // copy the information
            interp->Copy(gblock->dlvl(), gblock->shift(), &bsrc_neighbor, data_src, block_trg, data_trg);
        }
    }

    //................................................
    const bool do_coarse = (local_parent_.size() + ghost_parent_.size()) > 0;
    if (do_coarse) {
        // reset the coarse memory
        memset(coarse_ptr_, 0, interp->CoarseMemSize());
        // we use the full neighbor span
        const SubBlock bsrc_neighbor(M_GS, M_STRIDE, 0, M_N);
        //................................................
        // RMA the sibligns ghosts to the tmp
        for (const auto gblock : ghost_sibling_) {
            MPI_Aint disp_src  = gblock->data_src();
            rank_t   disp_rank = gblock->rank();
            // get the associated coarse block
            SubBlock block_trg;
            interp->CoarseFromFine(gblock, &block_trg);
            data_ptr data_trg = coarse_ptr_ + m_zeroidx(0, &block_trg);
            // interpolate, the level is 1 coarser and the shift is unchanged
            m_assert((gblock->dlvl() + 1) == 1, "the difference of level MUST be 1 or 0");
            interp->GetRma((gblock->dlvl() + 1), gblock->shift(), &bsrc_neighbor, disp_src, &block_trg, data_trg, disp_rank, mirrors_window);
        }
        // RMA the parent ghosts to the tmp
        for (const auto gblock : ghost_parent_) {
            MPI_Aint disp_src  = gblock->data_src();
            rank_t   disp_rank = gblock->rank();
            // get the associated coarse block
            SubBlock block_trg;
            interp->CoarseFromFine(gblock, &block_trg);
            data_ptr data_trg = coarse_ptr_ + m_zeroidx(0, &block_trg);
            // interpolate, the level is 1 coarser and the shift is unchanged
            m_assert((gblock->dlvl() + 1) == 0, "the difference of level MUST be 1 or 0");
            interp->GetRma((gblock->dlvl() + 1), gblock->shift(), &bsrc_neighbor, disp_src, &block_trg, data_trg, disp_rank, mirrors_window);
        }
        //................................................
        // copy the siblings to coarse
        for (const auto gblock : local_sibling_) {
            GridBlock* ngh_block = gblock->data_src();
            SubBlock   block_trg;
            interp->CoarseFromFine(gblock, &block_trg);
            data_ptr data_src = ngh_block->data(field, ida);
            data_ptr data_trg = coarse_ptr_ + m_zeroidx(0, &block_trg);
            // interpolate, the level is 1 coarser and the shift is unchanged
            m_assert((gblock->dlvl() + 1) == 1, "the difference of level MUST be 1");
            interp->Copy(gblock->dlvl() + 1, gblock->shift(), &bsrc_neighbor, data_src, &block_trg, data_trg);
        }
        // copy the parents to coarse
        for (const auto gblock : local_parent_) {
            GridBlock* ngh_block = gblock->data_src();
            SubBlock   block_trg;
            interp->CoarseFromFine(gblock, &block_trg);
            data_ptr data_src = ngh_block->data(field, ida);
            data_ptr data_trg = coarse_ptr_ + m_zeroidx(0, &block_trg);
            // interpolate, the level is 1 coarser and the shift is unchanged
            m_assert((gblock->dlvl() + 1) == 0, "the difference of level MUST be 1");
            interp->Copy(gblock->dlvl() + 1, gblock->shift(), &bsrc_neighbor, data_src, &block_trg, data_trg);
        }
    }
    //-------------------------------------------------------------------------
}

/**
 * @brief Do the second part of the GhostGet
 * 
 * - wait for the com to be over
 * - add the physical BC to my coarse temp
 * - interpolate to my ghost points
 * 
 * @param field 
 * @param ida 
 * @param interp 
 */
void GridBlock::GhostGet_Wait(const Field* field, const lda_t ida, const Wavelet* interp) {
    //-------------------------------------------------------------------------
    const bool do_coarse = (local_parent_.size() + ghost_parent_.size()) > 0;
    if (do_coarse) {
        //................................................
        {
            const SubBlock coarse_block(interp->CoarseNGhostFront(), interp->CoarseStride(), 0, M_HN);
            // copy myself to the coarse, one point out of 2
            const lid_t    shift[3] = {0, 0, 0};
            const data_ptr data_src = data(field, ida);
            const data_ptr data_trg = coarse_ptr_ + m_zeroidx(0, &coarse_block);
            // interpolate
            interp->Copy(1, shift, this, data_src, &coarse_block, data_trg);
        }

        //................................................
        // do here some physics, to completely fill the coarse block before the interpolation
        for (auto gblock : phys_) {
            // get the direction and the corresponding bctype
            const bctype_t bctype = field->bctype(ida, gblock->iface());
            // in the face direction, the start and the end are already correct, only the fstart changes
            SubBlock coarse_block;
            interp->CoarseFromFine(gblock, &coarse_block);
            lid_t fstart[3];
            interp->CoarseFromFine(face_start[gblock->iface()], fstart);
            data_ptr data_trg = coarse_ptr_ + m_zeroidx(0, &coarse_block);
            // apply the BC
            if (bctype == M_BC_NEU) {
                NeumanBoundary<M_WAVELET_N - 1> bc;
                bc(gblock->iface(), fstart, hgrid_, 0.0, &coarse_block, data_trg);
            } else if (bctype == M_BC_DIR) {
                DirichletBoundary<M_WAVELET_N - 1> bc;
                bc(gblock->iface(), fstart, hgrid_, 0.0, &coarse_block, data_trg);
            } else if (bctype == M_BC_EXTRAP) {
                ExtrapBoundary<M_WAVELET_N> bc;
                bc(gblock->iface(), fstart, hgrid_, 0.0, &coarse_block, data_trg);
            } else if (bctype == M_BC_ZERO) {
                ZeroBoundary bc;
                bc(gblock->iface(), fstart, hgrid_, 0.0, &coarse_block, data_trg);
            } else {
                m_assert(false, "this type of BC is not implemented yet or not valid %d", bctype);
            }
        }
        //................................................
        // refine from the coarse to the parents
        {
            // take the full coarse block and set the info in my GP
            const SubBlock block_src(interp->CoarseNGhostFront(), interp->CoarseStride(), -interp->CoarseNGhostFront(), M_HN + interp->CoarseNGhostBack());
            const data_ptr data_src = coarse_ptr_ + m_zeroidx(0, &block_src);
            lid_t          shift[3] = {0, 0, 0};

            for (const auto gblock : local_parent_) {
                interp->Interpolate(-1, shift, &block_src, data_src, gblock, data(field, ida));
            }
            for (const auto gblock : ghost_parent_) {
                interp->Interpolate(-1, shift, &block_src, data_src, gblock, data(field, ida));
            }
        }
    }
    //-------------------------------------------------------------------------
}

/**
 * @brief Do the first part of GhostPut: start the RMA put communications
 * 
 * - apply the physical BC to the best of my knowledge
 * - coarsen my data to obtain the scaling coefficient (interpolate, not copy!)
 * - put the obtained values to my neighbor's ghost points
 * 
 * @param field 
 * @param ida 
 * @param interp 
 * @param mirrors_window 
 */
void GridBlock::GhostPut_Post(const Field* field, const lda_t ida, const Wavelet* interp, MPI_Win mirrors_window) {
    //-------------------------------------------------------------------------
    const bool do_coarse = (local_parent_.size() + ghost_parent_.size()) > 0;
    if (do_coarse) {
        //................................................
        // reset the tmp to use for the put operations
        memset(coarse_ptr_, 0, interp->CoarseMemSize());
        //................................................
        // apply the physics to the best of my knowledge
        // some ghosts are missing but this is okay
        data_ptr data_trg = data(field, ida);
        for (auto gblock : phys_) {
            bctype_t bctype = field->bctype(ida, gblock->iface());
            // get the correct face_start
            if (bctype == M_BC_NEU) {
                NeumanBoundary<M_WAVELET_N - 1> bc;
                bc(gblock->iface(), face_start[gblock->iface()], hgrid_, 0.0, gblock, data_trg);
            } else if (bctype == M_BC_DIR) {
                DirichletBoundary<M_WAVELET_N - 1> bc;
                bc(gblock->iface(), face_start[gblock->iface()], hgrid_, 0.0, gblock, data_trg);
            } else if (bctype == M_BC_EXTRAP) {
                ExtrapBoundary<M_WAVELET_N> bc;
                bc(gblock->iface(), face_start[gblock->iface()], hgrid_, 0.0, gblock, data_trg);
            } else if (bctype == M_BC_ZERO) {
                ZeroBoundary bc;
                bc(gblock->iface(), face_start[gblock->iface()], hgrid_, 0.0, gblock, data_trg);
            } else {
                m_assert(false, "this type of BC is not implemented yet or not valid %d", bctype);
            }
        }
        {
            // //................................................
            // // I am now complete (except children GP), get my coarse representation
            const SubBlock coarse_block(interp->CoarseNGhostFront(), interp->CoarseStride(), 0, M_HN);
            data_ptr       data_coarse = coarse_ptr_ + m_zeroidx(0, &coarse_block);

            // the source block is the ghost extended block
            const lid_t    shift[3] = {0, 0, 0};
            const SubBlock me_extended(M_GS, M_STRIDE, -interp->nghost_front(), M_N + interp->nghost_back());
            // interpolate, the level is 1 coarser and the shift is unchanged
            interp->Interpolate(1, shift, &me_extended, data(field, ida), &coarse_block, data_coarse);

            //................................................
            // start the put commands to set the value to my neighbors
            // loop on the ghost list
            for (const auto gblock : ghost_parent_reverse_) {
                MPI_Aint disp_trg = gblock->data_src();
                rank_t   trg_rank = gblock->rank();
                // interpolate, the parent's mirror have been created to act on the tmp
                m_assert(gblock->dlvl() == 0, "we must have a level 0, here %d", gblock->dlvl());
                interp->PutRma(gblock->dlvl(), gblock->shift(), &coarse_block, coarse_ptr_, gblock, disp_trg, trg_rank, mirrors_window);
            }
            for (const auto gblock : local_parent_reverse_) {
                GridBlock* ngh_block = gblock->data_src();
                data_ptr   data_trg  = ngh_block->data(field, ida);
                // interpolate, the level is 1 coarser and the shift is unchanged
                m_assert(gblock->dlvl() == 0, "we must have a level 0, here %d", gblock->dlvl());
                interp->Copy(gblock->dlvl(), gblock->shift(), &coarse_block, data_coarse, gblock, data_trg);
            }
        }
    }
    //-------------------------------------------------------------------------
}

/**
 * @brief Do the second part of GhostPut
 * 
 * - the comm are over, I now have everything I need to compute the physical BC
 * 
 * @param field 
 * @param ida 
 * @param interp 
 */
void GridBlock::GhostPut_Wait(const Field* field, const lda_t ida, const Wavelet* interp) {
    //-------------------------------------------------------------------------
    data_ptr data_trg = data(field, ida);
    for (auto gblock : phys_) {
        bctype_t bctype = field->bctype(ida, gblock->iface());
        // get the correct face_start
        if (bctype == M_BC_NEU) {
            NeumanBoundary<M_WAVELET_N - 1> bc;
            bc(gblock->iface(), face_start[gblock->iface()], hgrid_, 0.0, gblock, data_trg);
        } else if (bctype == M_BC_DIR) {
            DirichletBoundary<M_WAVELET_N - 1> bc;
            bc(gblock->iface(), face_start[gblock->iface()], hgrid_, 0.0, gblock, data_trg);
        } else if (bctype == M_BC_EXTRAP) {
            ExtrapBoundary<M_WAVELET_N> bc;
            bc(gblock->iface(), face_start[gblock->iface()], hgrid_, 0.0, gblock, data_trg);
        } else if (bctype == M_BC_ZERO) {
            ZeroBoundary bc;
            bc(gblock->iface(), face_start[gblock->iface()], hgrid_, 0.0, gblock, data_trg);
        } else {
            m_assert(false, "this type of BC is not implemented yet or not valid %d", bctype);
        }
        //-------------------------------------------------------------------------
    }
    //-------------------------------------------------------------------------
}

// /**
//  * @brief Downsample the block and re-evaluate the boundary conditions on the coarse version
//  * 
//  * @param field 
//  * @param ida 
//  * @param interp 
//  */
// void GridBlock::Coarse_DownSampleWithBoundary(const Field* field, const lda_t ida, const Wavelet* interp, SubBlock* coarse_block) {
//     //-------------------------------------------------------------------------
//     // reset the tmp value
//     memset(coarse_ptr_, 0, interp->CoarseStride());
//     // get the coarse and extended SubBlocks
//     const lid_t n_ghost_front = (interp->ncriterion_front() + 1) / 2;
//     const lid_t n_ghost_back  = (interp->ncriterion_back() + 1) / 2;  // if we have 3 ghosts, we need 2 coarse and not 1...
//     const lid_t stride        = n_ghost_back + n_ghost_front + M_HN;
//     coarse_block->Reset(n_ghost_front, stride, -n_ghost_front, M_HN + n_ghost_back);
//     m_assert(stride <= interp->CoarseStride(), "ohoh the Coarse block is too small: %d <= %ld", stride, interp->CoarseStride());

//     const SubBlock extended_src(M_GS, M_STRIDE, -M_GS, M_N + M_GS);
//     const lid_t    shift[3] = {0, 0, 0};
//     const data_ptr data_trg = coarse_ptr_ + m_zeroidx(0, coarse_block);
//     const data_ptr data_src = data(field, ida);

//     // interpolate
//     interp->Copy(1, shift, &extended_src, data_src, coarse_block, data_trg);

//     //................................................
//     // over-write the BC
//     for (auto gblock : phys_) {
//         m_assert(false, "we shouldn't be here");
//         // get the direction and the corresponding bctype
//         const bctype_t bctype = field->bctype(ida, gblock->iface());
//         // in the face direction, the start and the end are already correct, only the fstart changes
//         SubBlock coarse_block;
//         interp->CoarseFromFine(gblock, &coarse_block);
//         lid_t fstart[3];
//         interp->CoarseFromFine(face_start[gblock->iface()], fstart);
//         // apply the BC
//         if (bctype == M_BC_NEU) {
//             NeumanBoundary<M_WAVELET_N - 1> bc;
//             bc(gblock->iface(), fstart, hgrid_, 0.0, &coarse_block, data_trg);
//         } else if (bctype == M_BC_DIR) {
//             DirichletBoundary<M_WAVELET_N - 1> bc;
//             bc(gblock->iface(), fstart, hgrid_, 0.0, &coarse_block, data_trg);
//         } else if (bctype == M_BC_EXTRAP) {
//             ExtrapBoundary<M_WAVELET_N> bc;
//             bc(gblock->iface(), fstart, hgrid_, 0.0, &coarse_block, data_trg);
//         } else if (bctype == M_BC_ZERO) {
//             ZeroBoundary bc;
//             bc(gblock->iface(), fstart, hgrid_, 0.0, &coarse_block, data_trg);
//         } else {
//             m_assert(false, "this type of BC is not implemented yet or not valid %d", bctype);
//         }
//     }
//     //-------------------------------------------------------------------------
// }