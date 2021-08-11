#ifndef SRC_TOOLS_TOOLSMPI_HPP_
#define SRC_TOOLS_TOOLSMPI_HPP_

#include "core/memdata.hpp"
#include "core/memspan.hpp"
#include "core/types.hpp"

void ToMPIDatatype(const MemLayout& layout, const MemSpan& span, const bidx_t scale, MPI_Datatype* xyz_type);

void GetRma(const level_t dlvl, const bidx_t shift[3],
            const MemLayout& layout_src, const MemSpan& span_src, MPI_Aint disp_src, rank_t src_rank,
            const MemLayout& layout_trg, const MemSpan& span_trg, MemData& data_trg, MPI_Win win);

void PutRma(const level_t dlvl, const bidx_t shift[3],
            const MemLayout& layout_src, const MemSpan& span_src, ConstMemData& data_src,
            const MemLayout& layout_trg, const MemSpan& span_trg, MPI_Aint disp_trg, rank_t trg_rank, MPI_Win win);

#endif