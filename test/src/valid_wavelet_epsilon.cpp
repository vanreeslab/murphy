#include <mpi.h>

#include "blas.hpp"
#include "doop.hpp"
#include "error.hpp"
#include "grid.hpp"
#include "gridcallback.hpp"
#include "gtest/gtest.h"
#include "ioh5.hpp"
#include "murphy.hpp"
#include "setvalues.hpp"
#include "subblock.hpp"
#include "wavelet.hpp"

#define DOUBLE_TOL 1e-13
#define BLVL 1

using std::list;

class valid_Wavelet_Epsilon : public ::testing::Test {
   protected:
    void SetUp() override{};
    void TearDown() override{};
};

static constexpr lid_t m_gs     = 8;
static constexpr lid_t m_n      = M_N;
static constexpr lid_t m_hn     = M_HN;
static constexpr lid_t m_stride = m_n + 6 * m_gs;

static InterpolatingWavelet* GetWavelet(const sid_t n, const sid_t nt) {
    if (n == 2 && nt == 2) {
        return (new Wavelet<2, 2>);
    } else if (n == 4 && nt == 0) {
        return (new Wavelet<4, 0>);
    } else if (n == 4 && nt == 2) {
        return (new Wavelet<4, 2>);
    } else if (n == 4 && nt == 4) {
        return (new Wavelet<4, 4>);
    } else if (n == 6 && nt == 0) {
        return (new Wavelet<6, 0>);
    } else if (n == 6 && nt == 2) {
        return (new Wavelet<6, 2>);
    } else if (n == 6 && nt == 4) {
        return (new Wavelet<6, 4>);
    }
    return nullptr;
};

//==============================================================================================================================
TEST_F(valid_Wavelet_Epsilon, epsilon_forced) {
    constexpr real_t hcoarse = 1.0 / (m_hn);
    constexpr real_t hfine   = 1.0 / (m_n);

    InterpolatingWavelet* interp = GetWavelet(M_WAVELET_N, M_WAVELET_NT);

    real_p ptr_fine   = (real_t*)m_calloc(m_stride * m_stride * m_stride * sizeof(real_t));
    real_p ptr_coarse = (real_t*)m_calloc(m_stride * m_stride * m_stride * sizeof(real_t));
    real_p ptr_tmp    = (real_t*)m_calloc(m_stride * m_stride * m_stride * sizeof(real_t));

    lid_t n[7][2] = {{2, 2}, {4, 0}, {4, 2}, {4, 4}, {6, 0}, {6, 2}, {6, 4}};

    // for (sid_t id = 0; id < 7; ++id) {
    // create the memory, start in full configuration
    SubBlock block_coarse(m_gs, m_hn + 2 * m_gs, -m_gs, m_hn + m_gs);
    SubBlock block_fine(3 * m_gs, m_n + 6 * m_gs, -3 * m_gs, m_n + 3 * m_gs);

    data_ptr data_fine   = ptr_fine + m_zeroidx(0, &block_fine);
    data_ptr data_coarse = ptr_coarse + m_zeroidx(0, &block_coarse);

    // fill the fine block!
    for (lid_t i2 = block_fine.start(2); i2 < block_fine.end(2); i2++) {
        for (lid_t i1 = block_fine.start(1); i1 < block_fine.end(1); i1++) {
            for (lid_t i0 = block_fine.start(0); i0 < block_fine.end(0); i0++) {
                real_t x      = i0 * hfine;
                real_t y      = i1 * hfine;
                real_t z      = i2 * hfine;
                real_t pos[3] = {x, y, z};

                data_fine[m_midx(i0, i1, i2, 0, &block_fine)] = cos(2 * M_PI * (x)) + cos(2 * M_PI * (y)) + cos(2 * M_PI * (z));
            }
        }
    }
    //................................................
    // compute the max detail
    block_fine.Reset(block_fine.gs(), block_fine.stride(), 0, m_n);
    // get the coarse version of life
    const lid_t m_hgs = m_gs / 2;
    SubBlock tmp_block(3 * m_hgs, m_hn + 6 * m_hgs, -3 * m_hgs, m_hn + 3 * m_hgs);
    data_ptr data_tmp = ptr_tmp + m_zeroidx(0, &tmp_block);
    for (lid_t i2 = tmp_block.start(2); i2 < tmp_block.end(2); i2++) {
        for (lid_t i1 = tmp_block.start(1); i1 < tmp_block.end(1); i1++) {
            for (lid_t i0 = tmp_block.start(0); i0 < tmp_block.end(0); i0++) {
                data_tmp[m_midx(i0, i1, i2, 0, &tmp_block)] = data_fine[m_midx(i0 * 2, i1 * 2, i2 * 2, 0, &block_fine)];
            }
        }
    }

    // InterpolatingWavelet* interp = GetWavelet(n[id][0], n[id][1]);
    real_t detail_max;
    interp->Details(&block_fine, data_fine, &tmp_block, data_tmp, &detail_max);

    //................................................
    // set the index for coarsening computation
    block_fine.Reset(block_fine.gs(), block_fine.stride(), -3 * m_gs, m_n + 3 * m_gs);
    // do the coarsening
    lid_t shift[3] = {0};
    interp->Interpolate(1, shift, &block_fine, data_fine, &block_coarse, data_coarse);

    // set the index for refinement computation
    block_fine.Reset(block_fine.gs(), block_fine.stride(), 0, m_n);
    interp->Interpolate(-1, shift, &block_coarse, data_coarse, &block_fine, data_fine);

    // now compute the error
    real_t erri = 0.0;
    for (lid_t i2 = block_fine.start(2); i2 < block_fine.end(2); i2++) {
        for (lid_t i1 = block_fine.start(1); i1 < block_fine.end(1); i1++) {
            for (lid_t i0 = block_fine.start(0); i0 < block_fine.end(0); i0++) {
                real_t x      = i0 * hfine;
                real_t y      = i1 * hfine;
                real_t z      = i2 * hfine;
                real_t pos[3] = {x, y, z};

                real_t exact = cos(2 * M_PI * (x)) + cos(2 * M_PI * (y)) + cos(2 * M_PI * (z));
                real_t value = data_fine[m_midx(i0, i1, i2, 0, &block_fine)];

                real_t err = fabs(exact - value);

                // if (err > erri) {
                //     m_log("update the value, I have %e instead of %e @ %d %d %d", value, exact, i0, i1, i2);
                // }
                erri = m_max(erri, err);
            }
        }
    }
    ASSERT_LE(erri, detail_max);
    m_log("[%s] erri = %e vs detail = %e", interp->Identity().c_str(), erri, detail_max);
    delete (interp);
    // }

    // ciao
    m_free(ptr_fine);
    m_free(ptr_coarse);

    ASSERT_NEAR(0.0, 0.0, DOUBLE_TOL);
}

//==============================================================================================================================
TEST_F(valid_Wavelet_Epsilon, epsilon_periodic_test) {
    // adapt the mesh
    real_t epsilon[2] = {1e-2, 1e-3};
    // real_t epsilon[1] = {1.0e-1};
    // lda_t  ieps       = 0;
    for (lda_t ieps = 0; ieps < 1; ++ieps) {
        level_t max_level   = 4;
        bool    periodic[3] = {true, true, true};
        lid_t   L[3]        = {1, 1, 1};
        Grid    grid(max_level, periodic, L, MPI_COMM_WORLD, nullptr);

        // create the field + the solution
        Field vort("vort", 3);
        grid.AddField(&vort);

        // create a patch of the current domain
        list<Patch> patch;
        real_t      origin[3] = {0.0, 0.0, 0.0};
        real_t      length[3] = {L[0], L[1], L[2]};
        patch.push_back(Patch(origin, length, max_level));

        const real_t         center[3] = {L[0] / 2.0, L[1] / 2.0, L[2] / 2.0};
        const lda_t          normal    = 2;
        const real_t         sigma     = 0.1;
        const real_t         radius    = 0.3;
        const real_t         cutoff    = (center[0] - radius) * 0.9;
        SetCompactVortexRing vr_init(normal, center, sigma, radius, cutoff);
        SetCompactVortexRing vr_init_full(normal, center, sigma, radius, cutoff, grid.interp());

        vr_init(&grid, &vort);

        grid.SetTol(1e+5, epsilon[ieps]);

        // do the coarsening, go the the min level if needed
        for (level_t il = max_level; il > 2; --il) {
            grid.Coarsen(&vort);
        }

        // go up again by forcing the refinement based on the patch
        for (level_t sil = 2; sil < max_level; ++sil) {
            grid.GhostPull(&vort);
            grid.Adapt(reinterpret_cast<void*>(&patch), nullptr, nullptr, &cback_Patch, &cback_Interpolate);
        }
        m_log("\t re-adaptation done! we have block between %d and %d", grid.MinLevel(), grid.MaxLevel());

        // recreate the solution
        Field sol("sol", 3);
        grid.AddField(&sol);
        vr_init(&grid, &sol);
        // compute the error
        real_t          err2, erri;
        ErrorCalculator error;
        error.Norms(&grid, &vort, &sol, &err2, &erri);
        m_log("==> error after reconstruction: epsilon %e: err2 = %e, erri = %e", epsilon[ieps], err2, erri);

        grid.DeleteField(&sol);

        ASSERT_LE(err2, epsilon[ieps]);
        ASSERT_LE(erri, epsilon[ieps]);
    }
}

//==============================================================================================================================
TEST_F(valid_Wavelet_Epsilon, epsilon_extrap_test) {
    // adapt the mesh
    real_t epsilon[3] = {1e-2, 1e-3};
    // real_t epsilon[1] = {1.0e-1};
    // lda_t  ieps       = 0;
    for (lda_t ieps = 0; ieps < 2; ++ieps) {
        level_t max_level   = 4;
        bool    periodic[3] = {false, false, false};
        lid_t   L[3]        = {1, 1, 1};
        Grid    grid(max_level, periodic, L, MPI_COMM_WORLD, nullptr);

        // create the field + the solution
        Field vort("vort", 3);
        grid.AddField(&vort);

        vort.bctype(M_BC_EXTRAP);

        // create a patch of the current domain
        list<Patch> patch;
        real_t      origin[3] = {0.0, 0.0, 0.0};
        real_t      length[3] = {L[0], L[1], L[2]};
        patch.push_back(Patch(origin, length, max_level));

        const real_t center[3] = {L[0] / 2.0, L[1] / 2.0, L[2] / 2.0};
        const lda_t  normal    = 2;
        const real_t sigma     = 0.05;
        const real_t radius    = 0.35;
        // const real_t  cutoff    = (center[0] - radius) * 0.9;
        SetVortexRing vr_init(normal, center, sigma, radius);
        SetVortexRing vr_init_full(normal, center, sigma, radius, grid.interp());

        vr_init(&grid, &vort);

        grid.SetTol(1e+5, epsilon[ieps]);

        // do the coarsening, go the the min level if needed
        for (level_t il = max_level; il > 2; --il) {
            grid.Coarsen(&vort);
        }

        grid.GhostPull(&vort);
        IOH5        io("data_test");
        std::string name = "vort_coarse" + std::to_string(epsilon[ieps]);
        io(&grid, &vort, name);

        // go up again by forcing the refinement based on the patch
        for (level_t sil = 2; sil < max_level; ++sil) {
            grid.GhostPull(&vort);
            grid.Adapt(reinterpret_cast<void*>(&patch), nullptr, nullptr, &cback_Patch, &cback_Interpolate);

            // recreate the solution
            Field sol("sol", 3);
            grid.AddField(&sol);
            vr_init(&grid, &sol);
            // compute the error
            real_t          err2, erri;
            ErrorCalculator error;
            error.Norms(&grid, &vort, &sol, &err2, &erri);
            m_log("==> error after reconstruction: epsilon %e: err2 = %e, erri = %e", epsilon[ieps], err2, erri);

            ASSERT_LE(err2, epsilon[ieps]);
            ASSERT_LE(erri, epsilon[ieps]);

            grid.DeleteField(&sol);
        }
        m_log("\t re-adaptation done! we have block between %d and %d", grid.MinLevel(), grid.MaxLevel());

        // recreate the solution
        Field sol("sol", 3);
        grid.AddField(&sol);
        vr_init(&grid, &sol);
        // compute the error
        real_t          err2, erri;
        ErrorCalculator error;
        error.Norms(&grid, &vort, &sol, &err2, &erri);
        m_log("==> error after reconstruction: epsilon %e: err2 = %e, erri = %e", epsilon[ieps], err2, erri);

        // grid.GhostPull(&vort);
        // IOH5 io("data_test");
        // io(&grid, &vort, "vort_fine");

        ASSERT_LE(err2, epsilon[ieps]);
        ASSERT_LE(erri, epsilon[ieps]);

        grid.DeleteField(&sol);
    }
}
