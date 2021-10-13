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
BlockDataType TypeToEnum<CartBlock>();// { return M_CARTBLOCK; }
template <>
BlockDataType TypeToEnum<GridBlock>();// { return M_GRIDBLOCK; }
template <>
BlockDataType TypeToEnum<const CartBlock>();// { return M_CARTBLOCK; }
template <>
BlockDataType TypeToEnum<const GridBlock>();// { return M_GRIDBLOCK; }

/**
 * @brief Check if a block type can be passed as an argument to a function
 */
static bool IsCompatibleBlockType(const BlockDataType type_to_convert_to, const BlockDataType type_to_convert_from) {
    m_begin;
    //--------------------------------------------------------------------------
    return (type_to_convert_to <= type_to_convert_from);
    //--------------------------------------------------------------------------
    m_end;
}

size_t BlockPointerSize(const BlockDataType bdt);
size_t BlockPartitionSize(const BlockDataType bdt, const lda_t lda);
void   AllocateBlockForP4est(const BlockDataType bdt, const real_t xyz[3], qdrt_t* p4est_quad);

#endif
