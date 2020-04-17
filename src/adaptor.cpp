// #include "adaptor.hpp"
// #include "gridcallback.hpp"
// #include "p8est_extended.h"

// Adaptor::Adaptor(ForestGrid* grid) {
//     m_begin;
//     //-------------------------------------------------------------------------
//     grid_ = grid;
//     //-------------------------------------------------------------------------
//     m_end;
// }

// void Adaptor::Refine(Interpolator* interp){
//     m_begin;
//     //-------------------------------------------------------------------------
//     tmp_interp_ = interp;

//     // compute the ghost needed by the interpolation
//         for (auto fid = fields_.begin(); fid != fields_.end(); fid++) {
//             // set the working field before entering the callback
//             // working_callback_field_ = fid->second;
//             // get the ghosts needed by the interpolation
//             GhostPull(fid->second);
//         }
//         // delete the soon-to be outdated ghost and mesh
//         delete (ghost_);
//         ResetP4estGhostMesh();
//         // set the grid in the forest for the callback
//         forest_->user_pointer = (void*)this;
//         // do the p4est interpolation by callback
//         p8est_refine_ext(forest_, 0, P8EST_MAXLEVEL, cback_Yes, nullptr, cback_Interpolate);
//         // balance the partition
//         // TODO
//         // create a new ghost and mesh
//         SetupP4estGhostMesh();
//         ghost_ = new Ghost(this);

//     //-------------------------------------------------------------------------
//     m_end;

// }