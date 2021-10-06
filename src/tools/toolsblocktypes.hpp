#ifndef SRC_GRID_BLOCKTYPETOOLS_
#define SRC_GRID_BLOCKTYPETOOLS_

#include "core/macros.hpp"
#include "core/types.hpp"
#include "grid/cartblock.hpp"
#include "grid/gridblock.hpp"
#include "grid/iimblock.hpp"

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
template <>
BlockDataType TypeToEnum<IIMBlock>() { return M_IIMBLOCK; }

// template <>
// BlockDataType TypeToEnum<const CartBlock>() { return M_CARTBLOCK; }
// template <>
// BlockDataType TypeToEnum<const GridBlock>() { return M_GRIDBLOCK; }
// template <>
// BlockDataType TypeToEnum<const IIMBlock>() { return M_IIMBLOCK; }

/**
 * @brief Check if a block type can be passed as an argument to a function
 */
bool IsCompatibleBlockType(BlockDataType type_to_convert_to, BlockDataType type_to_convert_from) {
    return (type_to_convert_to <= type_to_convert_from);
}

/**
 * @brief provide size of block pointer based on BlockDataType enum
 */
size_t BlockPointerSize(BlockDataType bdt) {
    m_assert(bdt != M_NULLTYPE, "cannot provide size of unspecified block type.");
    if (bdt == M_CARTBLOCK) return sizeof(CartBlock*);
    if (bdt == M_GRIDBLOCK) return sizeof(GridBlock*);
    if (bdt == M_IIMBLOCK) return sizeof(IIMBlock*);

    // if you get here, there is an error
    m_assert(false, "BlockDataType value not recognized.");
    return -1;
}

// /**
//  * @brief allocate a block and return a GridBlock pointer to its location
//  */
// template <typename B>
// B* AllocateGridBlock(BlockDataType bdt, const real_t length, const real_t xyz[3], const sid_t level);

#endif
