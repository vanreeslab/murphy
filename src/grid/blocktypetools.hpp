#ifndef SRC_GRID_BLOCKTYPETOOLS_
#define SRC_GRID_BLOCKTYPETOOLS_

#include "grid/cartblock.hpp"
#include "grid/gridblock.hpp"
#include "grid/iimblock.hpp"

/**
 * @brief Convert a block type to its BlockDataType enum value 
 */
template<typename BlockType> BlockDataType TypeToEnum(); 

/**
 * @brief Check if a block type can be passed as an argument to a function
 */
bool IsCompatibleBlockType(BlockDataType arg_type, BlockDataType block_type);

/**
 * @brief provide size of block pointer based on BlockDataType enum
 */
size_t BlockPointerSize(BlockDataType bdt);

/**
 * @brief allocate a block and return a GridBlock pointer to its location
 */
GridBlock* AllocateGridBlock(BlockDataType bdt, const real_t length, const real_t xyz[3], const sid_t level);

#endif


