#include <cmath>
#include <list>

#include "core/macros.hpp"
#include "core/types.hpp"
#include "grid/field.hpp"
#include "grid/grid.hpp"
#include "gtest/gtest.h"
#include "operator/error.hpp"
#include "operator/xblas.hpp"
#include "time/rk3_tvd.hpp"
#include "time/rkfunctor.hpp"
#include "tools/ioh5.hpp"

#define DOUBLE_TOL 1e-13

//=====================================================================================================
// small wrapper class around a Setvalue to access it using RK
class RKRhs : public RKFunctor, public SetValue {
   protected:
    bool accumulate_ = false;

    lda_t  ida_       = 0;
    real_t center_[3] = {0.0, 0.0, 0.0};
    real_t sigma_[3]  = {0.0, 0.0, 0.0};
    real_t alpha_     = 1.0;
    real_t time_      = 0.0;

   public:
    RKRhs(const real_t center[3], const real_t sigma[3], const real_t alpha, const lda_t ida, m_ptr<const Wavelet> interp) : SetValue(interp) {
        m_begin;
        m_assert(0 <= ida && ida <= 2, "the direction must be in [0,3[");
        //-------------------------------------------------------------------------
        for (int id = 0; id < 3; id++) {
            center_[id] = center[id];
            sigma_[id]  = sigma[id];
        }
        alpha_ = alpha;
        ida_   = ida;
        //-------------------------------------------------------------------------
        m_end;
    }

    // rk functor
    real_t cfl_rk3() const override { return 1.0; }
    real_t rdiff() const override { return 1.0; }

    void RhsSet(m_ptr<const Grid> grid, const real_t time, m_ptr<Field> field_u, m_ptr<Field> field_y) override {
        accumulate_ = false;
        time_       = time;
        m_log("time is now %e, evaluate the rhs with ida = %d", time_, ida_);
        m_log("evaluation in field %s", field_y->name().c_str());
        SetValue::operator()(grid, field_y);
        accumulate_ = false;
    };
    void RhsAcc(m_ptr<const Grid> grid, const real_t time, m_ptr<Field> field_u, m_ptr<Field> field_y) override {
        accumulate_ = true;
        time_       = time;
        SetValue::operator()(grid, field_y);
        accumulate_ = false;
    };

   protected:
    void FillGridBlock(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, m_ptr<Field> fid) override {
        //-------------------------------------------------------------------------
        real_t        pos[3];
        const real_t* xyz   = block->xyz();
        const real_t* hgrid = block->hgrid();

        real_t sigma     = sqrt(sigma_[0] * sigma_[0] + sigma_[1] * sigma_[1] + sigma_[2] * sigma_[2]);
        real_t oo_sigma2 = (sigma > 0.0) ? 1.0 / (sigma * sigma) : 0.0;
        real_t fact      = alpha_ * sqrt(1.0 / M_PI * oo_sigma2) * (-2.0 / pow(sigma_[ida_], 2));  // * alpha_ * sqrt(1.0 / M_PI) * (1.0 / (sigma * sigma * sigma));

        real_t reset = (accumulate_) ? 1.0 : 0.0;

        real_t center[3] = {center_[0] + time_ * (0 == ida_),
                            center_[1] + time_ * (1 == ida_),
                            center_[2] + time_ * (2 == ida_)};

        data_ptr block_data = block->data(fid);
        for (lda_t ida = ida_start_; ida < ida_end_; ida++) {
            real_t* data = block_data.Write(ida, block());
            for (lid_t i2 = start_; i2 < end_; i2++) {
                for (lid_t i1 = start_; i1 < end_; i1++) {
                    for (lid_t i0 = start_; i0 < end_; i0++) {
                        // get the position
                        m_pos(pos, i0, i1, i2, hgrid, xyz);
                        // compute the gaussian
                        const real_t rhox = (sigma_[0] > 0) ? ((pos[0] - center[0]) / sigma_[0]) : 0.0;
                        const real_t rhoy = (sigma_[1] > 0) ? ((pos[1] - center[1]) / sigma_[1]) : 0.0;
                        const real_t rhoz = (sigma_[2] > 0) ? ((pos[2] - center[2]) / sigma_[2]) : 0.0;
                        const real_t rho  = rhox * rhox + rhoy * rhoy + rhoz * rhoz;

                        const real_t dpos = (pos[ida_] - center[ida_]);

                        data[m_idx(i0, i1, i2)] = 0.0;  // *= reset;
                        data[m_idx(i0, i1, i2)] -= fact * std::exp(-rho) * dpos;
                    }
                }
            }
        }
        //-------------------------------------------------------------------------
    }
};

class valid_RK : public ::testing::Test {
    void SetUp() override{};
    void TearDown() override{};
};

using std::list;
using std::string;

/**
 * @brief test the RK3 implementation on a regular, uniform grid with a analytical advection source term
 * 
 */
TEST_F(valid_RK, rk3_tvd) {
    // setup the mesh
    bool  period[3] = {true, true, true};
    lid_t L[3]      = {3, 3, 3};

    // for (level_t il = 1; il < 3; ++il) {
    // get a uniform grid
    level_t     il = 1;
    Grid        grid(il, period, L, MPI_COMM_WORLD, nullptr);
    real_t      origin2[3] = {0.0, 0.0, 0.0};
    real_t      length2[3] = {L[0], L[1], L[2]};
    Patch       p2(origin2, length2, il);
    list<Patch> patch{p2};
    grid.Adapt(&patch);

    // obtain the time at evaluation as M_N/2 mesh points
    const real_t cfl          = 0.5;  // cfl = dt/dx
    const real_t dt_ref       = cfl * grid.FinestH();
    const iter_t iter_max_ref = 1;
    const real_t time_sol     = iter_max_ref * dt_ref;

    real_t erri_tvd[2] = {0.0, 0.0};

    // first we do iter_max, then iter_max*2
    for (short_t id = 0; id < 2; ++id) {
        const real_t dt       = dt_ref / (id + 1);
        const iter_t iter_max = iter_max_ref * (id + 1);

        m_log("test number %d, dt = %e, iter_max = %d", id, dt, iter_max);

        // get the fields
        string phiName = "phi" + std::to_string(id);
        string solName = "sol" + std::to_string(id);
        Field  phi(phiName, 1);
        Field  sol(solName, 1);
        grid.AddField(&phi);
        grid.AddField(&sol);

        // set the field - inital condition
        const real_t   center[3] = {1.5, 1.5, 1.5};
        const real_t   sigma[3]  = {0.1, 0.1, 0.1};
        const real_t   alpha     = 1.0;
        SetExponential phi_init(center, sigma, alpha, grid.interp());
        phi_init(&grid, &phi);

        // check the min and max
        real_t  min_init, max_init;
        BMinMax minmax;
        minmax(&grid, &phi, &min_init, &max_init);

        // set the rhs
        RKRhs rhs(center, sigma, alpha, 2, grid.interp());

        // set the RK
        RK3_TVD rk3(&grid, &phi, &rhs, nullptr);

        // integrate in time
        real_t time = 0.0;
        for (iter_t iter = 0; iter < iter_max; ++iter) {
            // do the dt
            rk3.DoDt(dt, &time);

            m_log("time is now %e", time);
            // check the dt value
            m_assert(time == (iter + 1) * dt, "the time must be now the current time-step");
        }

        m_log("solution is taken at time %e", time);
        // set the solution, assume velocity in Z: u = (0,0,1)
        const real_t center_sol[3] = {1.5, 1.5, 1.5 + 1 * time};
        // const real_t   center_sol[3] = {1.5, 1.5, 1.5};
        SetExponential sol_init(center_sol, sigma, alpha, grid.interp());
        sol_init(&grid, &sol);

        // compute the error
        ErrorCalculator error;
        error.Normi(&grid, &phi, &sol, erri_tvd + id);
        m_log("RK3 - TVD: checking iter_max = %d, ei = %e", iter_max, erri_tvd[id]);

        // check the min and max
        real_t min_final, max_final;
        minmax(&grid, &phi, &min_final, &max_final);

        // test the TVD
        m_log("init: %e to %e -> final %e to %e", min_init, max_init, min_final, max_final);
        // ASSERT_LE(min_init,min_final);
        // ASSERT_GE(max_init,max_final);

        // destroy the world
        grid.DeleteField(&phi);
        grid.DeleteField(&sol);
    }

    real_t convi = -log(erri_tvd[1] / erri_tvd[0]) / log(2);
    m_log("RK3 - TVD: the convergence orders are: norm_i:%e", convi);
    ASSERT_GE(convi, 3 - 0.1);
}