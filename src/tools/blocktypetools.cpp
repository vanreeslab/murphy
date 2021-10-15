// #include "tools/blocktypetools.hpp"


// template<typename BlockType> BlockDataType TypeToEnum() {
//     m_assert(false, "Block type not supported.");
//     return M_NULLTYPE;
// };

// template<> BlockDataType TypeToEnum<CartBlock>() { return M_CARTBLOCK; }
// template<> BlockDataType TypeToEnum<GridBlock>() { return M_GRIDBLOCK; }
// template<> BlockDataType TypeToEnum<IIMBlock>()  { return M_IIMBLOCK;  }

// template<> BlockDataType TypeToEnum<CartBlock const>() { return M_CARTBLOCK; }
// template<> BlockDataType TypeToEnum<GridBlock const>() { return M_GRIDBLOCK; }
// template<> BlockDataType TypeToEnum<IIMBlock  const>() { return M_IIMBLOCK;  }

// bool IsCompatibleBlockType(BlockDataType arg_type, BlockDataType block_type) {
//     m_assert((arg_type != M_NULLTYPE) && (block_type != M_NULLTYPE), "null block type not supported.");
//     return arg_type >= 
//     // if (arg_type == M_CARTBLOCK) {
//     //     return true;
//     // if (arg_type == M_GRIDBLOCK && block_type != M_CARTBLOCK)
//     //     return true;
//     // } else {
//     //     return arg_type == block_type;
//     // }
// }

// size_t BlockPointerSize(BlockDataType bdt) {
//     m_assert(bdt != M_NULLTYPE, "cannot provide size of unspecified block type.");
//     if (bdt == M_GRIDBLOCK) return sizeof(GridBlock*);
//     if (bdt == M_IIMBLOCK)  return sizeof(IIMBlock*);

//     // if you get here, there is an error
//     m_assert(false, "BlockDataType value not recognized.");
//     return -1; 
// }

// GridBlock* AllocateGridBlock(BlockDataType bdt, const real_t length, const real_t xyz[3], const sid_t level) {
//     m_assert(bdt != M_NULLTYPE, "cannot allocate unspecified block type.");
//     m_assert(bdt != M_CARTBLOCK, "cannot pass a CartBlock via GridBlock*.");
//     if (bdt == M_GRIDBLOCK) return (GridBlock*)(new GridBlock(length, xyz, level));
//     if (bdt == M_IIMBLOCK)  return (GridBlock*)(new IIMBlock(length, xyz, level));

//     // if you get here, there is an error
//     m_assert(false, "BlockDataType value not recognized.");
//     return nullptr; 
// }

