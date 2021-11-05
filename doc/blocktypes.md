# Managing block types in MURPHY
In MURPHY, the type of block managed by a grid can be specified at runtime. To make this possible, 
each block type is given a label of enumerated type `BlockDataType`. Operations that create blocks, destroy blocks, or cast between block pointer types use this enum to track type information and perform type checking.

The downside to this approach is that adding a new block type to MURPHY requires, depending on the case, modifying both the `BlockDataType` enumeration and all functions that make use of it. The full process is as follows:
 - Create a header defining a new block type that inherits from `CartBlock` or one of its child classes. If the block will require any kind of ghosting operation, it should inherit from `GridBlock` or one of its child classes.
 - Override the methods `CartBlock::PartitionDataOffset()`, `CartBlock::PartitionDataPack()`, and `CartBlock::PartitionDataUnPack()`. This will allow the partitioner to transfer extra block information between ranks.
 - Add a label for the new block type to the `BlockDataType` enumeration in `core/types.hpp`.
 - Include the new block header in `tools/blocktypetools.hpp`. Then add the new block type/enum to all of the template specializations in `tools/blocktypetools.hpp`and all of the implementations in `tools/toolsblocktypes.cpp`.
 - Add the new block type/enum to the implementation of `get_cback_CreateBlock()` and `get_cback_DestroyBlock` in `grid/gridcallback.cpp`.
 - If needed, add a new template specialization `Stencil<~new block type~>` to `operator/stencil.cpp`.
