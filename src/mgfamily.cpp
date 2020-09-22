#include "mgfamily.hpp"
#include "subblock.hpp"
#include <p8est_bits.h>

MGFamily::MGFamily(const lid_t num_children) {
    m_begin;
    //-------------------------------------------------------------------------
    parent_count_ = 0;
    m_assert(num_children % P8EST_CHILDREN == 0, "the number of children = %d must be a multiple of %d", num_children, P8EST_CHILDREN);
    parents_  = reinterpret_cast<GridBlock**>(m_calloc(sizeof(GridBlock*) * (num_children / P8EST_CHILDREN)));
    children_ = reinterpret_cast<GridBlock**>(m_calloc(sizeof(GridBlock*) * num_children));
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief Destroy the MGFamily and the parents as nobody knows that they have been created
 * 
 */
MGFamily::~MGFamily() {
    for (lid_t id = 0; id < parent_count_; id++) {
        delete (parents_[id]);
    }
    // I need to free the parents I created as nobody else know that I have created them
    m_free(parents_);
    m_free(children_);
}

void MGFamily::AddMembers(GridBlock* parent, GridBlock* children[P8EST_CHILDREN]) {
    //-------------------------------------------------------------------------
    parents_[parent_count_] = parent;
    for (sid_t ic = 0; ic < P8EST_CHILDREN; ic++) {
        children_[parent_count_ * P8EST_CHILDREN + ic] = children[ic];
    }
    parent_count_++;
    //-------------------------------------------------------------------------
}

/**
 * @brief compute the coarsen of the children to the parents: field_trg = alpha * field_cst + interp(field_src)
 * 
 * @param field_src a fine field that needs to be interpolated
 * @param field_trg a coarse field that will contain the interpolation + alpha * field_cst
 * @param field_cst a constant coarse field to add to the interpolation
 * @param alpha the multiplication factor for the constant field
 * @param interp the interpolator to use
 */
void MGFamily::ToChildren(Field* field_src, Field* field_trg, Field* field_cst, const real_t alpha, Interpolator* interp) {
    m_begin;
    m_assert(field_src->lda() == field_trg->lda(), "the source and traget dimensions MUST match");
    //-------------------------------------------------------------------------
    // create an empty SubBlock representing the valid GP
    lid_t     parent_start[3] = {0, 0, 0};
    lid_t     parent_end[3]   = {0, 0, 0};
    SubBlock* mem_block       = new SubBlock(M_GS, M_STRIDE, parent_start, parent_end);

    for (lid_t id = 0; id < parent_count_; id++) {
        GridBlock* parent = parents_[id];
        for (sid_t ic = 0; ic < P8EST_CHILDREN; ic++) {
            GridBlock* child = children_[id * P8EST_CHILDREN + ic];

            // get the shift given the child position (trg) wrt to the parent framework (src)
            lid_t shift[3];
            shift[0] = (child->xyz(0) - parent->xyz(0)) / parent->hgrid(0);
            shift[1] = (child->xyz(1) - parent->xyz(1)) / parent->hgrid(1);
            shift[2] = (child->xyz(2) - parent->xyz(2)) / parent->hgrid(2);

            // create a subblock that contains the correct memory zone, we refine
            for (int id = 0; id < 3; id++) {
                parent_start[id] = shift[id] - M_GS;
                parent_end[id]   = shift[id] + M_HN + M_GS;
            }
            mem_block->Reset(M_GS, M_STRIDE, parent_start, parent_end);

            // interpolate for every dimension
            for (sid_t ida = 0; ida < field_src->lda(); ida++) {
                // get the pointers
                // interp->Interpolate(-1, shift, mem_block, parent->data(field_src, ida), child, child->data(field_trg, ida), alpha, child->data(field_cst, ida));
            }
        }
    }
    // un-validate the ghost status
    field_trg->ghost_status(false);
    // delete the temp subblock
    delete(mem_block);
    //-------------------------------------------------------------------------
    m_end;
}

void MGFamily::ToParents(Field* field_src, Field* field_trg, Field* field_cst, const real_t alpha, Interpolator* interp) {
    m_begin;
    m_assert(field_src->lda() == field_trg->lda(),"the source and traget dimensions MUST match");
    //-------------------------------------------------------------------------
    // create an empty SubBlock representing the valid GP
    lid_t     parent_start[3] = {0, 0, 0};
    lid_t     parent_end[3]   = {0, 0, 0};
    SubBlock* mem_block       = new SubBlock(M_GS, M_STRIDE, parent_start, parent_end);

    for (lid_t id = 0; id < parent_count_; id++) {
        GridBlock* parent = parents_[id];
        for (sid_t ic = 0; ic < P8EST_CHILDREN; ic++) {
            GridBlock* child = children_[id * P8EST_CHILDREN + ic];

            // get the shift given the parent position (trg) wrt to the child framework (src)
            lid_t shift[3];
            shift[0] = (parent->xyz(0) - child->xyz(0)) / child->hgrid(0);
            shift[1] = (parent->xyz(1) - child->xyz(1)) / child->hgrid(1);
            shift[2] = (parent->xyz(2) - child->xyz(2)) / child->hgrid(2);

            // create a subblock that contains the correct memory zone, we coarsen
            parent_start[0] = (child->xyz(0) - parent->xyz(0)) / parent->hgrid(0);
            parent_start[1] = (child->xyz(1) - parent->xyz(1)) / parent->hgrid(1);
            parent_start[2] = (child->xyz(2) - parent->xyz(2)) / parent->hgrid(2);
            for (int id = 0; id < 3; id++) {
                parent_end[id] = parent_start[id] + M_HN;
            }
            mem_block->Reset(M_GS, M_STRIDE, parent_start, parent_end);

            // interpolate for every dimension
            for (sid_t ida = 0; ida < field_src->lda(); ida++) {
                // get the pointers
                // interp->Interpolate(+1, shift, child, child->data(field_src, ida), mem_block, parent->data(field_trg, ida), alpha, parent->data(field_cst, ida));
            }
        }
    }
    // un-validate the ghost status
    field_trg->ghost_status(false);
    // delete the temp memblock
    delete(mem_block);
    //-------------------------------------------------------------------------
    m_end;
}
