#include <cmath>

#include "boundary.hpp"
#include "error.hpp"
#include "field.hpp"
#include "gtest/gtest.h"
#include "ioh5.hpp"
#include "laplacian.hpp"
#include "murphy.hpp"
#include "setvalues.hpp"
#include "subblock.hpp"
#include "wavelet.hpp"

#define DOUBLE_TOL 1e-9

class valid_Stencil : public ::testing::Test {
   protected:
    Grid* grid_;

    int  lvl_         = 1;
    bool periodic_[3] = {false, false, false};
    int  l_[3]        = {1, 1, 1};

    Field* vort_;
    Field* diff_;
    Field* sol_;

    MemPool* mem_pool_;

    void SetUp() override {
        mem_pool_ = new MemPool();
        grid_     = new Grid(lvl_, periodic_, l_, MPI_COMM_WORLD, nullptr, mem_pool_);

        //add the field
        vort_ = new Field("vorticity", 1);
        diff_ = new Field("diffusion", 1);
        sol_  = new Field("solution", 1);
        grid_->AddField(vort_);
        grid_->AddField(diff_);
        grid_->AddField(sol_);
    };
    void TearDown() override {
        delete (vort_);
        delete (diff_);
        delete (sol_);
        delete (grid_);
        delete(mem_pool_);
    };
};

TEST_F(valid_Stencil, laplacian_o2) {
    vort_->bctype(M_BC_EXTRAP_3);

    // x^2 + y^2 + z^2
    real_t     dir[3]  = {1.0, 1.0, 1.0};
    lid_t      deg[3]  = {2, 2, 2};
    SetPolynom polynom = SetPolynom(deg, dir);
    polynom(grid_, vort_);

    // 2 + 2 + 2
    real_t     dir2[3]  = {2.0, 2.0, 2.0};
    lid_t      deg2[3]  = {0, 0, 0};
    SetPolynom polynom2 = SetPolynom(deg2, dir2);
    polynom2(grid_, sol_);

    LaplacianCross<3> lapla = LaplacianCross<3>();
    lapla(grid_, vort_, diff_);

    real_t          norm2, normi;
    ErrorCalculator myerr;
    myerr.Norms(grid_, diff_, sol_, &norm2, &normi);

    // ASSERT_LE(norm2, normi);
    ASSERT_NEAR(norm2, 0.0, DOUBLE_TOL);
    ASSERT_NEAR(normi, 0.0, DOUBLE_TOL);
}

TEST_F(valid_Stencil, laplacian_o4) {
    vort_->bctype(M_BC_EXTRAP_5);

    // x^4 + y^4 + z^4
    real_t     dir[3]  = {1.0, 1.0, 1.0};
    lid_t      deg[3]  = {4, 4, 4};
    SetPolynom polynom = SetPolynom(deg, dir);
    polynom(grid_, vort_);

    // 12 x^2 + 12 y^2 + 12 z^2
    real_t     dir2[3]  = {12.0, 12.0, 12.0};
    lid_t      deg2[3]  = {2, 2, 2};
    SetPolynom polynom2 = SetPolynom(deg2, dir2);
    polynom2(grid_, sol_);

    LaplacianCross<5> lapla = LaplacianCross<5>();
    lapla(grid_, vort_, diff_);

    real_t          norm2, normi;
    ErrorCalculator myerr;
    myerr.Norms(grid_, diff_, sol_, &norm2, &normi);

    // ASSERT_LE(norm2, normi);
    ASSERT_NEAR(norm2, 0.0, DOUBLE_TOL);
    ASSERT_NEAR(normi, 0.0, DOUBLE_TOL);
}

#define N_CONV 2

TEST_F(valid_Stencil, convergence_laplacian_o2_boundary3) {
    vort_->bctype(M_BC_EXTRAP_3);

    real_t erri[N_CONV] = {0.0};
    real_t origin1[3]   = {0.0, 0.0, 0.0};
    real_t length1[3]   = {1.0 * l_[0], 1.0 * l_[1], 0.5 * l_[2]};
    real_t origin2[3]   = {0.0, 0.0, 0.5 * l_[2]};
    real_t length2[3]   = {1.0 * l_[0], 1.0 * l_[1], 0.5 * l_[2]};

    real_t normi[N_CONV];

    for (int il = lvl_; il < lvl_ + N_CONV; il++) {
        list<Patch> patches;
        patches.push_back(Patch(origin1, length1, il));
        patches.push_back(Patch(origin2, length2, il + 1));
        grid_->Adapt(&patches);

        // x^3 + y^4 + z^3
        real_t     dir[3]  = {1.0, 1.0, 1.0};
        lid_t      deg[3]  = {3, 4, 3};
        SetPolynom polynom = SetPolynom(deg, dir);
        polynom(grid_, vort_);

        // 2 + 2 + 2
        real_t     dir2[3]  = {6.0, 12.0, 6.0};
        lid_t      deg2[3]  = {1, 2, 1};
        SetPolynom polynom2 = SetPolynom(deg2, dir2);
        polynom2(grid_, sol_);

        LaplacianCross<3> lapla = LaplacianCross<3>();
        lapla(grid_, vort_, diff_);

        real_t          norm2;
        ErrorCalculator myerr;
        myerr.Norms(grid_, diff_, sol_, &norm2, normi + (il - lvl_));
        if (grid_->mpirank() == 0) {
            printf("lvl = %d error = %e %e\n", il, norm2, normi[il - lvl_]);
        }
    }
    if (grid_->mpirank() == 0) {
        printf("tested convergence order: %e vs %e\n", normi[N_CONV - 2] / normi[N_CONV - 1], 1.95);
    }
    ASSERT_GT(normi[N_CONV - 2] / normi[N_CONV - 1], 1.95);
}

TEST_F(valid_Stencil, convergence_laplacian_o2_boundary4) {
    vort_->bctype(M_BC_EXTRAP_4);

    real_t erri[N_CONV] = {0.0};
    real_t origin1[3]   = {0.0, 0.0, 0.0};
    real_t length1[3]   = {1.0 * l_[0], 1.0 * l_[1], 0.5 * l_[2]};
    real_t origin2[3]   = {0.0, 0.0, 0.5 * l_[2]};
    real_t length2[3]   = {1.0 * l_[0], 1.0 * l_[1], 0.5 * l_[2]};

    real_t normi[N_CONV];

    for (int il = lvl_; il < lvl_ + N_CONV; il++) {
        list<Patch> patches;
        patches.push_back(Patch(origin1, length1, il));
        patches.push_back(Patch(origin2, length2, il + 1));
        grid_->Adapt(&patches);

        // x^3 + y^4 + z^3
        real_t     dir[3]  = {1.0, 1.0, 1.0};
        lid_t      deg[3]  = {3, 4, 3};
        SetPolynom polynom = SetPolynom(deg, dir);
        polynom(grid_, vort_);

        // 2 + 2 + 2
        real_t     dir2[3]  = {6.0, 12.0, 6.0};
        lid_t      deg2[3]  = {1, 2, 1};
        SetPolynom polynom2 = SetPolynom(deg2, dir2);
        polynom2(grid_, sol_);

        LaplacianCross<3> lapla = LaplacianCross<3>();
        lapla(grid_, vort_, diff_);

        real_t          norm2;
        ErrorCalculator myerr;
        myerr.Norms(grid_, diff_, sol_, &norm2, normi + (il - lvl_));
        if (grid_->mpirank() == 0) {
            printf("lvl = %d error = %e %e\n", il, norm2, normi[il - lvl_]);
        }
    }
    if (grid_->mpirank() == 0) {
        printf("tested convergence order: %e vs %e\n", normi[N_CONV - 2] / normi[N_CONV - 1], 3.99);
    }
    ASSERT_GT(normi[N_CONV - 2] / normi[N_CONV - 1], 3.99);
}

TEST_F(valid_Stencil, convergence_laplacian_o4_boundary4) {
    vort_->bctype(M_BC_EXTRAP_4);

    real_t erri[N_CONV] = {0.0};
    real_t origin1[3]   = {0.0, 0.0, 0.0};
    real_t length1[3]   = {1.0 * l_[0], 1.0 * l_[1], 0.5 * l_[2]};
    real_t origin2[3]   = {0.0, 0.0, 0.5 * l_[2]};
    real_t length2[3]   = {1.0 * l_[0], 1.0 * l_[1], 0.5 * l_[2]};

    real_t normi[N_CONV];

    for (int il = lvl_; il < lvl_ + N_CONV; il++) {
        list<Patch> patches;
        patches.push_back(Patch(origin1, length1, il));
        patches.push_back(Patch(origin2, length2, il + 1));
        grid_->Adapt(&patches);

        // x^6 + y^5 + z^7
        real_t     dir[3]  = {1.0, 1.0, 1.0};
        lid_t      deg[3]  = {6, 5, 7};
        SetPolynom polynom = SetPolynom(deg, dir);
        polynom(grid_, vort_);

        // solution
        real_t     dir2[3]  = {30.0, 20.0, 42.0};
        lid_t      deg2[3]  = {4, 3, 5};
        SetPolynom polynom2 = SetPolynom(deg2, dir2);
        polynom2(grid_, sol_);

        LaplacianCross<5> lapla = LaplacianCross<5>();
        lapla(grid_, vort_, diff_);

        real_t          norm2;
        ErrorCalculator myerr;
        myerr.Norms(grid_, diff_, sol_, &norm2, normi + (il - lvl_));
        if (grid_->mpirank() == 0) {
            printf("lvl = %d error = %e %e\n", il, norm2, normi[il - lvl_]);
        }
    }
    if (grid_->mpirank() == 0) {
        printf("tested convergence order: %e vs %e\n", normi[N_CONV - 2] / normi[N_CONV - 1], 3.82);
    }
    ASSERT_GT(normi[N_CONV - 2] / normi[N_CONV - 1], 3.82);
}

TEST_F(valid_Stencil, convergence_laplacian_o4_boundary5) {
    vort_->bctype(M_BC_EXTRAP_5);

    real_t erri[N_CONV] = {0.0};
    real_t origin1[3]   = {0.0, 0.0, 0.0};
    real_t length1[3]   = {1.0 * l_[0], 1.0 * l_[1], 0.5 * l_[2]};
    real_t origin2[3]   = {0.0, 0.0, 0.5 * l_[2]};
    real_t length2[3]   = {1.0 * l_[0], 1.0 * l_[1], 0.5 * l_[2]};

    real_t normi[N_CONV];

    for (int il = lvl_; il < lvl_ + N_CONV; il++) {
        list<Patch> patches;
        patches.push_back(Patch(origin1, length1, il));
        patches.push_back(Patch(origin2, length2, il + 1));
        grid_->Adapt(&patches);

        // x^6 + y^5 + z^7
        real_t     dir[3]  = {1.0, 1.0, 1.0};
        lid_t      deg[3]  = {6, 5, 7};
        SetPolynom polynom = SetPolynom(deg, dir);
        polynom(grid_, vort_);

        // solution
        real_t     dir2[3]  = {30.0, 20.0, 42.0};
        lid_t      deg2[3]  = {4, 3, 5};
        SetPolynom polynom2 = SetPolynom(deg2, dir2);
        polynom2(grid_, sol_);

        LaplacianCross<5> lapla = LaplacianCross<5>();
        lapla(grid_, vort_, diff_);

        real_t          norm2;
        ErrorCalculator myerr;
        myerr.Norms(grid_, diff_, sol_, &norm2, normi + (il - lvl_));
        if (grid_->mpirank() == 0) {
            printf("lvl = %d error = %e %e\n", il, norm2, normi[il - lvl_]);
        }
    }
    if (grid_->mpirank() == 0) {
        printf("tested convergence order: %e vs %e\n", normi[N_CONV - 2] / normi[N_CONV - 1], 7.77);
    }
    ASSERT_GT(normi[N_CONV - 2] / normi[N_CONV - 1], 7.77);
}
