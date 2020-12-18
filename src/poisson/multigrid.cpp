// #include "multigrid.hpp"

// #include <mpi.h>

// #include <iostream>
// #include <limits>

// #include "daxpy.hpp"
// #include "core/macros.hpp"
#include "core/types.hpp"
// #include "error.hpp"
// #include "gaussseidel.hpp"
// #include "gridcallback.hpp"
// #include "ioh5.hpp"
// #include "jacobi.hpp"
// #include "laplacian.hpp"

// using std::numeric_limits;

// // #define MG_GAUSSSEIDEL
// #define MG_JACOBI_ORDER 3

// Multigrid::Multigrid(Grid* grid, sid_t fft_level, Field* rhs, Field* sol, Field* res) {
//     m_begin;
//     //-------------------------------------------------------------------------
//     // store the desired fft level and which field will be used as what
//     fft_level_ = fft_level;
//     m_assert(grid->IsAField(rhs), "the source field MUST exist in the grid");
//     fields_nickname_["rhs"]  = rhs->name();
//     map_fields_[rhs->name()] = rhs;
//     m_assert(grid->IsAField(sol), "the solution field MUST exist in the grid");
//     map_fields_[sol->name()] = sol;
//     fields_nickname_["sol"]  = sol->name();
//     m_assert(grid->IsAField(res), "the residual field MUST exist in the grid");
//     map_fields_[res->name()] = res;
//     fields_nickname_["res"]  = res->name();
//     m_assert(rhs->lda() == sol->lda(), "the dimension between the source and the solution MUST match");
//     m_assert(rhs->lda() == res->lda(), "the dimension between the source and the residual MUST match");

//     // get by how many level I have to do
//     sid_t l_max_level = p4est_MaxLocalLevel(grid->p4est_forest());
//     // upate the counters on every rank
//     m_assert(sizeof(l_max_level) == 1, "change the MPI datatype bellow");
//     MPI_Allreduce(&l_max_level, &max_level_, 1, MPI_CHAR, MPI_MAX, MPI_COMM_WORLD);
//     n_level_ = max_level_ - fft_level_;

//     // allocate the grids and the partitionner
//     grids_    = reinterpret_cast<Grid**>(m_calloc(sizeof(Grid*) * (n_level_ + 1)));
//     families_ = reinterpret_cast<MGFamily**>(m_calloc(sizeof(MGFamily*) * n_level_));
//     parts_    = reinterpret_cast<Partitioner**>(m_calloc(sizeof(Partitioner*) * n_level_));

//     // remember the original grid
//     grids_[n_level_] = grid;

//     // coarsen the grid
//     for (sid_t il = (n_level_ - 1); il >= 0; il--) {
//         // create a new grid and get the old information
//         grids_[il] = new Grid();
//         grids_[il]->CopyFrom(grids_[il + 1]);

//         // regester myself in the user_pointer
//         tmp_ilevel_               = il;
//         p8est_t* curr_forest      = grids_[il]->p4est_forest();
//         curr_forest->user_pointer = this;

//         // get the number of CHILDREN on the fine level, i.e. the number of children entering the family
//         sc_array_t quadarray = grids_[il + 1]->p4est_mesh()->quad_level[fft_level + il + 1];
//         m_assert(quadarray.elem_count < numeric_limits<lid_t>::max(), "the number of quad is too big");
//         // create the new family
//         families_[il] = new MGFamily(quadarray.elem_count);

//         // coarsen by one level the quads on the highest level and only allocate the fields in the MG map
//         p8est_coarsen_ext(curr_forest, 0, 0, cback_Level, nullptr, cback_MGCreateFamilly);
//         m_assert(families_[il]->parent_count() == (quadarray.elem_count / P8EST_CHILDREN), "those two numbers must match: %d vs %ld", families_[il]->parent_count(), (quadarray.elem_count / P8EST_CHILDREN));

//         // partition the grid and remember it, only the fields in the multigrid map exist now!
//         parts_[il] = new Partitioner(&map_fields_, grids_[il], false);
//         // create a consistent Ghost layout
//         grids_[il]->SetupGhost();
//         // the grid is now partitioned on a coarser level with new ghosts
//     }

//     direct_solver_ = new FFTSolver(grids_[0], sol, fft_level_);

//     //-------------------------------------------------------------------------
//     m_end;
// }

// Multigrid::~Multigrid() {
//     m_begin;
//     //-------------------------------------------------------------------------

//     for (sid_t il = (n_level_ - 1); il >= 0; il--) {
//         delete (parts_[il]);
//         delete (families_[il]);
//         delete (grids_[il]);
//     }

//     m_free(parts_);
//     m_free(families_);
//     m_free(grids_);

//     delete direct_solver_;

//     //-------------------------------------------------------------------------
//     m_end;
// }

// void Multigrid::Solve() {
//     m_begin;
//     //-------------------------------------------------------------------------
//     // get the fields
//     Field* sol = map_fields_.at(fields_nickname_.at("sol"));
//     Field* res = map_fields_.at(fields_nickname_.at("res"));
//     Field* rhs = map_fields_.at(fields_nickname_.at("rhs"));
//     // get which field will be send from one level to another
//     map<string, Field*> map_sol_field;
//     map_sol_field[fields_nickname_.at("sol")] = sol;
//     map<string, Field*> map_rhs_field;
//     map_rhs_field[fields_nickname_.at("rhs")] = rhs;

//     // needed tools
//     Daxpy daxpy_minus = Daxpy(-1.0);
//     Daxpy daxpy_plus  = Daxpy(+1.0);

// #ifndef MG_GAUSSSEIDEL
//     Jacobi<MG_JACOBI_ORDER>         jacobi = Jacobi<MG_JACOBI_ORDER>();
//     LaplacianCross<MG_JACOBI_ORDER> lapla  = LaplacianCross<MG_JACOBI_ORDER>();
// #else
//     LaplacianCross<3> lapla = LaplacianCross<3>();
//     GaussSeidel<3>    gs    = GaussSeidel<3>(1.0);
// #endif
// #ifndef NDEBUG
//     // needed to track the error
//     real_t          norm2;
//     ErrorCalculator error2 = ErrorCalculator();
// #endif

//     // string fname;
//     // IOH5   dump = IOH5("data");
//     // initial error
// #ifndef NDEBUG
//     lapla(grids_[n_level_], sol, res);
//     real_t normi;
//     error2.Norms(grids_[n_level_], res, rhs, &norm2, &normi);
//     m_log("initial error = %e and %e", norm2, normi);
// #endif
//     for (sid_t icycle = 0; icycle < 5; icycle++) {
//         m_log("---- cycle %d",icycle);
//         //---------------------
//         // GOING DOWN
//         for (sid_t il = n_level_; il > 0; il--) {
//             // this is the current grid
//             Grid* grid = grids_[il];
//             // Smooth that stuff
//             for (sid_t ie = 0; ie < eta_1_; ie++) {
// #ifndef MG_GAUSSSEIDEL
//                 jacobi(sol, rhs, res, grid);
// #else
//                 gs(sol, rhs, nullptr, grid);
// #endif
//                 // #ifndef NDEBUG
//                 //         lapla(grid, sol, res);
//                 //         error2.Norms(grid, res, rhs, &norm2,&normi);
//                 //         m_log("down: **Jacobi %d, error = %e %e**", ie, norm2, normi);
//                 // #endif
//             }
//             // fname = "down_sol_level" + std::to_string(il);
//             // dump(grid, sol, fname);
//             // compute the residual as rhs - A x to send to the coarser level
//             // this operation is overkilled since it executed on the WHOLE active leaf, while only the finest level is concerned
//             lapla(grid, sol, res);
// #ifndef NDEBUG
//             error2.Norms(grid, res, rhs, &norm2, &normi);
//             m_log("down: leaving level %d, error = %e %e", il, norm2, normi);
// #endif
//             // fine residual = - fine residula + fine rhs
//             // this operation is overkilled since it executed on the WHOLE active leaf, while only the finest level is concerned
//             daxpy_minus(grid, res, rhs, res);

//             // do the interpolation from the residual to the rhs of the next level
//             // for the finer level: the coarser rhs = 0.0 * coarse rhs + interp(residual), the rest of the grid remain unchanged
//             families_[il - 1]->ToParents(res, rhs, rhs, 0.0, grid->interp());
//             // do the partitioning, send the new rhs
//             parts_[il - 1]->Start(&map_rhs_field, M_FORWARD);
//             parts_[il - 1]->End(&map_rhs_field, M_FORWARD);
//         }
//         //---------------------
//         // do the direct solve
//         Grid* grid_fft = grids_[0];
//         //     // do the jacobi becuse I don't have anything else for the moment
//         //     for (lid_t ie = 0; ie < 45; ie++) {
//         // #ifndef MG_GAUSSSEIDEL
//         //         jacobi(sol, rhs, res, grid_fft);
//         // #else
//         //         gs(sol, rhs, nullptr, grid_fft);
//         // #endif
//         //         // compute the residual as rhs - A x to send to the coarser level
//         //     }
//         //     IOH5 dump = IOH5("data");
//         //     dump(grid_fft, sol, "jacob_sol");
//         //     dump(grid_fft, rhs, "jacob_rhs");

//         // // do the FFT instead
//         (*direct_solver_)(grid_fft, rhs, sol);
//         // dump(grid_fft, sol, "fft_sol");
//         // dump(grid_fft, rhs, "fft_rhs");

// #ifndef NDEBUG
//         lapla(grid_fft, sol, res);
//         // dump(grid_fft, res, "fft_lapla");
//         error2.Norms(grid_fft, res, rhs, &norm2, &normi);
//         m_log("after the direct: error = %e %e", norm2, normi);
// #endif

//         //---------------------

//         // GOING UP
//         for (sid_t il = 1; il <= n_level_; il++) {
//             // get the ghost for the solution as we need them to refine
//             grids_[il - 1]->GhostPull(sol);
//             // fname = "sol_coarse_level_" + std::to_string(il-1);
//             // dump(grid, sol,fname);
//             // do the partitioning, send the solution + ghost
//             parts_[il - 1]->Start(&map_sol_field, M_BACKWARD);
//             parts_[il - 1]->End(&map_sol_field, M_BACKWARD);
//             // do the interpolation using the same interpolation as the one used to go down
//             // fname = "before_familly_level_" + std::to_string(il);
//             // dump(grids_[il-1], sol,fname);
//             // the finer solution = the existing solution + Interp( the coarse solution)
//             // families_[il - 1]->ToChildren(sol, res, res,0.0, grids_[il]->interp());
//             families_[il - 1]->ToChildren(sol, sol, sol, 1.0, grids_[il]->interp());
//             // fname = "after_familly_level_" + std::to_string(il);
//             // dump(grids_[il],sol,fname);
//             // we are now on the new grid!!
//             Grid* grid = grids_[il];
//             // update the solution with it's coarsegrid version
//             // daxpy_plus(grid, res, sol, sol);  // sol = res + sol
//             // fname = "sol_coarse_interp_level_" + std::to_string(il);
//             // dump(grid, res,fname);
//             // fname = "sol_level_" + std::to_string(il);
//             // dump(grid, sol,fname);
// #ifndef NDEBUG
//             lapla(grid, sol, res);
//             // fname = "sol_level_lapla_" + std::to_string(il);
//             // dump(grid, res, fname);
//             error2.Norm2(grid, res, rhs, &norm2);
//             m_log("up: entering level %d, error = %e", il, norm2);
// #endif
//             // Gauss-Seidel - laplacian(sol) = rhs
//             for (sid_t ie = 0; ie < eta_2_; ie++) {
// #ifndef MG_GAUSSSEIDEL
//                 jacobi(sol, rhs, res, grid);
// #else
//                 gs(sol, rhs, nullptr, grid);
// #endif
//                 // #ifndef NDEBUG
//                 //             lapla(grid, sol, res);
//                 //             error2.Norms(grid, res, rhs, &norm2, &normi);
//                 //             m_log("up: Jacobi %d, error = `%16.16e %16.16e`", ie, norm2, normi);
//                 // #endif
//             }
// #ifndef NDEBUG
//             lapla(grid, sol, res);
//             error2.Norms(grid, res, rhs, &norm2, &normi);
//             m_log("up: leaving level %d, error = %e %e", il, norm2, normi);
// #endif
//         }
//     }
//     //-------------------------------------------------------------------------
//     m_end;
// }