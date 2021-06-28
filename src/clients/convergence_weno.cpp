#include "clients/convergence_weno.hpp"

#include <list>
#include <string>

#include "grid/field.hpp"
#include "grid/grid.hpp"
#include "operator/advection.hpp"
#include "operator/error.hpp"

using std::list;
using std::string;
using std::to_string;

#define N_PT 3

ConvergenceWeno::~ConvergenceWeno() {
    //-------------------------------------------------------------------------
    //-------------------------------------------------------------------------
}

void ConvergenceWeno::InitParam(ParserArguments* param) {
    //-------------------------------------------------------------------------
    ilevel_ = m_max(param->init_lvl, 1);
    adapt_  = !param->no_adapt;
    //-------------------------------------------------------------------------
}

void ConvergenceWeno::Run() {
    m_begin;
    //-------------------------------------------------------------------------
    // setup the mesh
    const bool  period[3] = {true, true, true};
    const lid_t L[3]      = {3, 3, 3};
    const real_t freq[3]    = {2.0, 2.0, 2.0};

    const real_t rand_vel[3] = {-1.0 + ((real_t)std::rand() / (real_t)RAND_MAX) * 2.0,
                                -1.0 + ((real_t)std::rand() / (real_t)RAND_MAX) * 2.0,
                                -1.0 + ((real_t)std::rand() / (real_t)RAND_MAX) * 2.0};

    for (level_t il = 0; il < N_PT; ++il) {
        Grid grid(il + ilevel_, period, L, MPI_COMM_WORLD, nullptr);

        // create the patch refinement to refine the middle tree
        real_t origin11[3] = {1.0, 1.0, 1.0};
        real_t length1[3]  = {1.0, 1.0, 1.0};
        Patch  p1(origin11, length1, il + ilevel_ + 1);
        // real_t origin2[3] = {0.0, 0.0, 0.0};
        // real_t length2[3] = {3.0, 3.0, 3.0};
        // Patch  p2(origin2, length2, il + ilevel_ - 1);

        // list<Patch> patch{p11,p2};
        list<Patch> patch{p1};
        if (adapt_) {
            grid.Adapt(&patch);
        }

        // create the needed fields
        string fieldName = "field" + std::to_string(il);
        // string solName   = "sol" + std::to_string(il);
        string velName   = "vel" + std::to_string(il);
        string diffName  = "deriv" + std::to_string(il);
        Field  test(fieldName, 1);
        // Field  sol(solName, 1);
        Field  vel(velName, 3);
        Field  dtest(diffName, 1);
        grid.AddField(&test);
        // grid.AddField(&sol);
        grid.AddField(&vel);
        grid.AddField(&dtest);

        // set a constant velocity
        {
            lambda_setvalue_t lambda_vel = [=](const bidx_t i0, const bidx_t i1, const bidx_t i2, const CartBlock* const block, const Field* const fid)->void{
                block->data(fid,0).Write(i0,i1,i2)[0] = rand_vel[0];
                block->data(fid,1).Write(i0,i1,i2)[0] = rand_vel[1];
                block->data(fid,2).Write(i0,i1,i2)[0] = rand_vel[2];
            };
            SetValue vel_init(lambda_vel,grid.interp());
            vel_init(&grid, &vel);
            // grid.GhostPull(&vel);
        }

        // set the test field
        {
            lambda_setvalue_t lambda_cos = [=](const bidx_t i0, const bidx_t i1, const bidx_t i2, const CartBlock* const block, const Field* const fid) -> void {
                real_t pos[3];
                block->pos(i0, i1, i2, pos);
                real_t* data = block->data(fid, 0).Write(i0, i1, i2);

                data[0] = std::cos(2.0 * M_PI * freq[0] / L[0] * pos[0]) +
                          std::cos(2.0 * M_PI * freq[0] / L[0] * pos[1]) +
                          std::cos(2.0 * M_PI * freq[0] / L[0] * pos[2]);
            };
            SetValue field_init(lambda_cos, grid.interp());
            field_init(&grid, &test);
        }

        // // set the solution
        // {
        //     // -> the solution: - u* df/dx - v * df/dy - w*df/dz
        //     const real_t freq[3] = {2.0, 2.0, 2.0};

        //     lambda_setvalue_t lambda_sin = [=](const bidx_t i0, const bidx_t i1, const bidx_t i2, const CartBlock* const block, const Field* const fid) -> void {
        //         real_t pos[3];
        //         block->pos(i0, i1, i2, pos);
        //         real_t* data = block->data(fid, 0).Write(i0, i1, i2);

        //         data[0] = (rand_vel[0] * 2.0 * M_PI * freq[0] / L[0]) * std::sin(2.0 * M_PI * freq[0] / L[0] * pos[0]) +
        //                   (rand_vel[1] * 2.0 * M_PI * freq[1] / L[1]) * std::sin(2.0 * M_PI * freq[0] / L[0] * pos[1]) +
        //                   (rand_vel[2] * 2.0 * M_PI * freq[2] / L[2]) * std::sin(2.0 * M_PI * freq[0] / L[0] * pos[2]);
        //     };
        //     SetValue sol_init(lambda_sin, grid.interp());
        //     sol_init(&grid, &sol);
        // }
        lambda_error_t lambda_sol = [=](const bidx_t i0, const bidx_t i1, const bidx_t i2, const CartBlock* const block) -> real_t {
            real_t pos[3];
            block->pos(i0, i1, i2, pos);

            return (rand_vel[0] * 2.0 * M_PI * freq[0] / L[0]) * std::sin(2.0 * M_PI * freq[0] / L[0] * pos[0]) +
                   (rand_vel[1] * 2.0 * M_PI * freq[1] / L[1]) * std::sin(2.0 * M_PI * freq[0] / L[0] * pos[1]) +
                   (rand_vel[2] * 2.0 * M_PI * freq[2] / L[2]) * std::sin(2.0 * M_PI * freq[0] / L[0] * pos[2]);
        };

        const bool do_weno_3 = grid.NGhostFront() >= 2 && grid.NGhostBack() >= 2;
        const bool do_weno_5 = grid.NGhostFront() >= 3 && grid.NGhostBack() >= 3;

        rank_t rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        real_t hmax = grid.FinestH();

        if (do_weno_3) {
            Advection<M_WENO_Z, 3> adv(&vel);
            adv(&grid, &test, &dtest);
            m_log("error weno 3");
            // now, we need to check
            
            Error  error;
            real_t erri, err2;
            error.Norms(&grid, &dtest, &lambda_sol, &err2, &erri);
            if (rank == 0) {
                string fname     = "data/conv_weno3_a" + to_string(adapt_) + "_w" + to_string(M_WAVELET_N) + to_string(M_WAVELET_NT) + ".data";
                FILE*  file_diag = fopen(fname.c_str(), "a+");
                fprintf(file_diag, "%e %e %e\n", hmax, err2, erri);
                fclose(file_diag);
            }
            m_log("%e %e %e", hmax, err2, erri);
        }
        if (do_weno_5) {
            Advection<M_WENO_Z, 5> adv(&vel);
            adv(&grid, &test, &dtest);
            m_log("error weno 3");
            // now, we need to check
            Error  error;
            real_t erri, err2;
            error.Norms(&grid, &dtest,  &lambda_sol, &err2, &erri);
            if (rank == 0) {
                string fname     = "data/conv_weno5_a" + to_string(adapt_) + "_w" + to_string(M_WAVELET_N) + to_string(M_WAVELET_NT) + ".data";
                FILE*  file_diag = fopen(fname.c_str(), "a+");
                fprintf(file_diag, "%e %e %e\n", hmax, err2, erri);
                fclose(file_diag);
            }
            m_log("%e %e %e", hmax, err2, erri);
        }

        // grid.DeleteField(&sol);
        grid.DeleteField(&vel);
        grid.DeleteField(&test);
        grid.DeleteField(&dtest);
    }
    //-------------------------------------------------------------------------
    m_end;
}
