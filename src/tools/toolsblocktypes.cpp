#include "toolsblocktypes.hpp"

template <>
BlockDataType TypeToEnum<CartBlock>() { return M_CARTBLOCK; }
template <>
BlockDataType TypeToEnum<GridBlock>() { return M_GRIDBLOCK; }
template <>
BlockDataType TypeToEnum<const CartBlock>() { return M_CARTBLOCK; }
template <>
BlockDataType TypeToEnum<const GridBlock>() { return M_GRIDBLOCK; }

/**
 * @brief provide size of block pointer based on BlockDataType enum
 */
size_t BlockPointerSize(const BlockDataType bdt) {
    m_begin;
    //--------------------------------------------------------------------------
    m_assert(bdt != M_NULLTYPE, "cannot provide size of unspecified block type.");
    if (bdt == M_CARTBLOCK) return sizeof(CartBlock*);
    if (bdt == M_GRIDBLOCK) return sizeof(GridBlock*);

    // if you get here, there is an error
    m_assert(false, "BlockDataType value not recognized.");
    return -1;
    //--------------------------------------------------------------------------
    m_end;
}

/**
 * @brief Return the size needed by one block during the partitioning
 * 
 * @warning we create a random block and query the size, it's not the best but it works fine enough
 */
size_t BlockPartitionSize(const BlockDataType bdt, const lda_t lda) {
    m_begin;
    //--------------------------------------------------------------------------
    m_assert(bdt != M_NULLTYPE, "cannot provide size of unspecified block type.");
    const real_t xyz[3] = {0.0, 0.0, 0.0};
    if (bdt == M_CARTBLOCK) {
        CartBlock    block(1.0, xyz, 0);
        const size_t n_elem   = block.BlockLayout().n_elem;
        const size_t n_offset = block.PartitionDataOffset();

        return n_elem * lda + n_offset;
    }
    if (bdt == M_GRIDBLOCK) {
        GridBlock    block(1.0, xyz, 0);
        const size_t n_elem   = block.BlockLayout().n_elem;
        const size_t n_offset = block.PartitionDataOffset();

        return n_elem * lda + n_offset;
    }

    // if you get here, there is an error
    m_assert(false, "BlockDataType value not recognized.");
    return 0;
    //--------------------------------------------------------------------------
    m_end;
}

/**
 * @brief allocate a block and return a GridBlock pointer to its location
 */
void AllocateBlockForP4est(const BlockDataType bdt, const real_t xyz[3], qdrt_t* p4est_quad) {
    //--------------------------------------------------------------------------
    const real_t length = p4est_QuadLen(p4est_quad->level);
    if (bdt == M_CARTBLOCK) {
        CartBlock* block = new CartBlock(length, xyz, p4est_quad->level);
        p4est_SetBlock(p4est_quad, block);
    }
    if (bdt == M_GRIDBLOCK) {
        GridBlock* block = new GridBlock(length, xyz, p4est_quad->level);
        p4est_SetBlock(p4est_quad, block);
    }
    //--------------------------------------------------------------------------
}