#ifndef SRC_GRID_BLOCKTYPETOOLS_
#define SRC_GRID_BLOCKTYPETOOLS_

#include "core/macros.hpp"
#include "core/types.hpp"
#include "grid/cartblock.hpp"
#include "grid/gridblock.hpp"
#include "grid/iimblock.hpp"
#include "toolsp4est.hpp"

/**
 * @brief Convert a block type to its BlockDataType enum value 
 */
template <typename BlockType>
BlockDataType TypeToEnum() {
    m_assert(false, "Block type not supported.");
    return M_NULLTYPE;
};

template <>
BlockDataType TypeToEnum<CartBlock>() { return M_CARTBLOCK; }
template <>
BlockDataType TypeToEnum<GridBlock>() { return M_GRIDBLOCK; }

/**
 * @brief Check if a block type can be passed as an argument to a function
 */
bool IsCompatibleBlockType(const BlockDataType type_to_convert_to, const BlockDataType type_to_convert_from) {
    m_begin;
    //--------------------------------------------------------------------------
    return (type_to_convert_to <= type_to_convert_from);
    //--------------------------------------------------------------------------
    m_end;
}

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
 * @brief allocate a block and return a GridBlock pointer to its location
 */
void AllocateBlockForP4est(const BlockDataType bdt, const real_t xyz[3], qdrt_t* p4est_quad) {
    //--------------------------------------------------------------------------
    const real_t length = p4est_QuadLen(p4est_quad->level);
    if (bdt == M_CARTBLOCK) {
        CartBlock* block = new CartBlock(length, xyz, p4est_quad->level);
        p4est_SetBlock(p4est_quad, block);
    }
    if (bdt == M_CARTBLOCK) {
        CartBlock* block = new CartBlock(length, xyz, p4est_quad->level);
        p4est_SetBlock(p4est_quad, block);
    }
    //--------------------------------------------------------------------------
}

#endif
