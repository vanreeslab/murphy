// #include "fftsolver.hpp"


// FLUPS_BoundaryType flups_get_bc(bool isper, bctype_t mybc)
// {
//     if (isper)
//     {
//         return PER;
//     }
//     else
//     {
//         switch (mybc)
//         {
//         case M_BC_EVEN:
//             return EVEN;
//             break;
//         case M_BC_ODD:
//             return ODD;
//             break;
//         case M_BC_ZERO:
//             return UNB;
//             break;
//         case M_BC_EXTRAP_3:
//             return UNB;
//             break;
//         case M_BC_EXTRAP_4:
//             return UNB;
//             break;
//         case M_BC_EXTRAP_5:
//             return UNB;
//             break;
//         default:
//             m_assert(false, "error");
//             break;
//         }
//     }
//     return NONE;
// };


// FFTSolver::FFTSolver(Grid *grid, Field *sol, const int fft_level) {
//     m_begin;
//     //-------------------------------------------------------------------------
//     // init the FLUPS Solver
//     const lid_t nglob[3] = {M_N * grid->domain_length(0), M_N * grid->domain_length(1), M_N * grid->domain_length(2)};
//     const lid_t nproc[3] = {grid->domain_length(0), grid->domain_length(1), grid->domain_length(2)};
//     m_assert(fft_level == 0, "for the moment we must use the root as the fft level");
//     flups_topo_ = flups_topo_new(0, sol->lda(), nglob, nproc, false, nullptr, M_ALIGNMENT, MPI_COMM_WORLD);

//     // get the needed infor
//     const real_t h[3] = {1.0 / M_N, 1.0 / M_N, 1.0 / M_N};
//     const real_t L[3] = {grid->domain_length(0), grid->domain_length(1), grid->domain_length(2)};

//     m_log("hfactor = %f %f %f",h[0],h[1],h[2]);
//     m_log("Lfactor = %f %f %f",L[0],L[1],L[2]);
//     m_log("nglob = %d %d %d",nglob[0],nglob[1],nglob[2]);
//     m_log("nproc = %d %d %d",nproc[0],nproc[1],nproc[2]);

//     for (int id = 0; id < 3; id++) {
//         for (int is = 0; is < 2; is++) {
//             flups_bc_[id][is] = (FLUPS_BoundaryType *)m_calloc(sizeof(int) * sol->lda());
//         }
//     }
//     m_log("is grid per? %d %d %d",grid->domain_periodic(0),grid->domain_periodic(1),grid->domain_periodic(2));
//     for (sid_t ida = 0; ida < sol->lda(); ida++) {
//         flups_bc_[0][0][ida] = flups_get_bc(grid->domain_periodic(0), sol->bctype(ida, 0));
//         flups_bc_[0][1][ida] = flups_get_bc(grid->domain_periodic(0), sol->bctype(ida, 1));
//         flups_bc_[1][0][ida] = flups_get_bc(grid->domain_periodic(1), sol->bctype(ida, 2));
//         flups_bc_[1][1][ida] = flups_get_bc(grid->domain_periodic(1), sol->bctype(ida, 3));
//         flups_bc_[2][0][ida] = flups_get_bc(grid->domain_periodic(2), sol->bctype(ida, 4));
//         flups_bc_[2][1][ida] = flups_get_bc(grid->domain_periodic(2), sol->bctype(ida, 5));
//     }
//     flups_solver_ = flups_init(flups_topo_, flups_bc_, h, L, NOD);

//     flups_set_greenType(flups_solver_, LGF_2);
//     flups_setup(flups_solver_, true);


//     data_rhs_ = reinterpret_cast<real_t*>(m_calloc(M_N*M_N*M_N*sol->lda()*sizeof(real_t)));
//     data_sol_ = reinterpret_cast<real_t*>(m_calloc(M_N*M_N*M_N*sol->lda()*sizeof(real_t)));
//     //-------------------------------------------------------------------------
//     m_end;
// }

// FFTSolver::~FFTSolver() {
//     m_begin;
//     //-------------------------------------------------------------------------
//     // flups
//     flups_cleanup(flups_solver_);
//     flups_topo_free(flups_topo_);
//     for (int id = 0; id < 3; id++) {
//         for (int is = 0; is < 2; is++) {
//             m_free(flups_bc_[id][is]);
//         }
//     }

//     m_free(data_rhs_);
//     m_free(data_sol_);

//     //-------------------------------------------------------------------------
//     m_end;
// }

// void FFTSolver::ApplyOpF2F(const qid_t *qid, GridBlock *block, Field *fid_src, Field *fid_trg) {
//     //-------------------------------------------------------------------------
//     // do the copy the rhs only
//     for (sid_t ida = 0; ida < fid_src->lda(); ida++) {
//         real_p fid_data = block->data(fid_src, ida);
//         m_assume_aligned(fid_data);
//         m_assume_aligned(data_rhs_);
//         for (lid_t i2 = 0; i2 < M_N; i2++) {
//             for (lid_t i1 = 0; i1 < M_N; i1++) {
//                 for (lid_t i0 = 0; i0 < M_N; i0++) {
//                     data_rhs_[i0 + M_N * (i1 + M_N * (i2 + M_N * ida))] = fid_data[m_idx(i0, i1, i2)];
//                 }
//             }
//         }
//     }

//     m_log("solving with flups");
//     // solve
//     flups_solve(flups_solver_, data_sol_, data_rhs_,STD);

//     //copy back the solution
//     for (sid_t ida = 0; ida < fid_src->lda(); ida++) {
//         real_p fid_data = block->data(fid_trg, ida);
//         m_assume_aligned(fid_data);
//         m_assume_aligned(data_rhs_);
//         for (lid_t i2 = 0; i2 < M_N; i2++) {
//             for (lid_t i1 = 0; i1 < M_N; i1++) {
//                 for (lid_t i0 = 0; i0 < M_N; i0++) {
//                     fid_data[m_idx(i0, i1, i2)] = data_sol_[i0 + M_N * (i1 + M_N * (i2 + M_N * ida))];
//                 }
//             }
//         }
//     }
//     //-------------------------------------------------------------------------
// }