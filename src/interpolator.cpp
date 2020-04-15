#include "interpolator.hpp"

// void Interpolator::operator()(const Field* field, const sid_t ida, GhostBlock* block){
//     // get the relation
//     const sid_t dlvl = block->dlvl();
//     const lid_t* shift = block->shift();

//     // get the trg info
//     real_p trg_data = block->origin()->data(field,ida);

//     if()
//     // create a subblock describing the ghost memory
//     int gstart[3] = {0,0,0};
//     int gend[3] = {M_N,M_N,M_N};
//     SubBlock ghostblock = SubBlock(0,M_N,gstart,gend);
//     real_p src_data = block->data_src();

//     interpolate_(dlvl,shift,&ghostblock,src_data,block,trg_data);
// }

// void Interpolator::operator()(const Field* field, const sid_t ida, GridBlock* block_trg, GridBlock* block_src){
//     // get the relation
//     const sid_t dlvl = 0;
//     const lid_t shift[3] ={0,0,0};

//     // get the trg and memory info
//     real_p src_data = block_src->data(field,ida);
//     real_p trg_data = block_trg->data(field,ida);

//     // interpolate!
//     interpolate_(dlvl,shift,block_src,src_data,block_trg,trg_data);
// }

// void Interpolator::operator()(const Field* field, const sid_t ida, GridBlock* block_trg, MemLayout* block_src, real_p data_src){
//     // get the relation
//     const sid_t dlvl = 0;
//     const lid_t shift[3] ={0,0,0};

//     // get the trg and memory info
//     real_p trg_data = block_trg->data(field,ida);

//     // interpolate!
//     interpolate_(dlvl,shift,block_src,data_src,block_trg,trg_data);
// }

/**
    * @brief interpolates the data from data_src to data_trg
    * 
    * @param dlvl the difference of level: level_src - level_trg, i.e. > 0 means coarsening, = 0 means copy and < 0 means refinement
    * @param shift the position of the trg (0,0,0) in the src framework (and resolution!)
    * @param block_src description of @ref data_src memory layout 
    * @param data_src the 0-position of the src memory, i.e. the memory location of (0,0,0) for the source
    * @param block_trg descripiton of @ref data_trg
    * @param data_trg the 0-position of the trg memory, i.e. the memory location of (0,0,0) for the target
    */
void Interpolator::Interpolate(const sid_t dlvl, const lid_t shift[3], MemLayout* block_src, real_p data_src, MemLayout* block_trg, real_p data_trg) {
    m_begin;
    m_assert(dlvl <= 2, "we cannot handle a difference in level > 2");
    m_assert(dlvl >= -1, "we cannot handle a level too coarse ");
    //-------------------------------------------------------------------------
    // create the interpolation context
    interp_ctx_t ctx;


    m_verb("entering interpolator with shift = %d %d %d",shift[0],shift[1],shift[2]);
    m_verb("entering interpolator with srcstart = %d %d %d",block_src->start(0),block_src->start(1),block_src->start(2));
    m_verb("entering interpolator with srcend = %d %d %d",block_src->end(0),block_src->end(1),block_src->end(2));
    m_verb("entering interpolator with trgstart = %d %d %d",block_trg->start(0),block_trg->start(1),block_trg->start(2));
    m_verb("entering interpolator with trgend = %d %d %d",block_trg->end(0),block_trg->end(1),block_trg->end(2));

    // get memory details
    for (int id = 0; id < 3; id++) {
        // the parent starting and ending is place form the child point of view
#ifndef NDEBUG
        // the starting position of the source contains the ghosts, by definition
        ctx.srcstart[id] = block_src->start(id) - shift[id];
        ctx.srcend[id]   = block_src->end(id) - shift[id];
#endif
        ctx.trgstart[id] = block_trg->start(id);
        ctx.trgend[id]   = block_trg->end(id);
    }
    ctx.srcstr = block_src->stride();
    ctx.trgstr = block_trg->stride();

    // get the correct aligned arrays
    // note: since the adresses refer to (0,0,0), we have a ghostsize of 0
    ctx.sdata = data_src + m_midx(shift[0], shift[1], shift[2], 0, block_src);
    ctx.tdata = data_trg;

    // call the correct function
    if (dlvl == -1) {
        Refine_(&ctx);
    } else if (dlvl == 0) {
        Copy_(&ctx);
    } else if (dlvl > 0) {
        Coarsen_(&ctx, dlvl);
    }
    //-------------------------------------------------------------------------
    m_end;
}