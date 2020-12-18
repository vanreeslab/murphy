// #ifndef SRC_FFTSOLVER_HPP_
// #define SRC_FFTSOLVER_HPP_

// #include "flups.h"
// #include "grid.hpp"
// #include "core/macros.hpp"
#include "core/types.hpp"
// #include "operator.hpp"

// class FFTSolver : public OperatorF2F {
//    protected:
//     // flups stuffs
//     FLUPS_Topology *    flups_topo_;
//     FLUPS_Solver *      flups_solver_;
//     FLUPS_BoundaryType *flups_bc_[3][2];


//     real_p data_rhs_ = nullptr;
//     real_p data_sol_ = nullptr;

//    public:
//     FFTSolver(Grid *grid, Field *sol, const int fft_level);
//     ~FFTSolver();

//     void ApplyOpF2F(const qid_t *qid, GridBlock *block, Field *fid_src, Field *fid_trg) override;
// };

// #endif  // SRC_FFTSOLVER_HPP_