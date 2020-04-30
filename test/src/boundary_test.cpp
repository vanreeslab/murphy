
#include "boundary.hpp"

#include <cmath>

#include "error.hpp"
#include "gtest/gtest.h"
#include "murphy.hpp"
#include "subblock.hpp"
#include "wavelet.hpp"

#define DOUBLE_TOL 1e-13

class valid_Boundary : public ::testing::Test {
   protected:
    SubBlock* block_;
    real_p    data_;
    real_t    hgrid_[3];

    lid_t start_[3];
    lid_t end_[3];

    PhysBlock* physblock_[6];

    void SetUp() override {
        data_ = (real_t*)m_calloc(M_STRIDE * M_STRIDE * M_STRIDE * sizeof(real_t));

        for (int id = 0; id < 3; id++) {
            start_[id] = 0;
            end_[id]   = M_N;
            hgrid_[id] = 1.0 / (M_N);
        }

        block_ = new SubBlock(M_GS, M_STRIDE, start_, end_);
        data_  = data_ + m_zeroidx(0, block_);
    };
    void TearDown() override {
        m_free(data_ - m_zeroidx(0, block_));
        delete (block_);
    };
};

//==============================================================================================================================
TEST_F(valid_Boundary, bc_even_left) {
    // fill the source
    for (int i2 = -M_GS; i2 < M_N + M_GS; i2++) {
        for (int i1 = -M_GS; i1 < M_N + M_GS; i1++) {
            for (int i0 = 0; i0 < M_N; i0++) {
                real_t x = ((real_t)i0 + 0.5) * hgrid_[0];
                real_t y = ((real_t)i1 + 0.5) * hgrid_[1];
                real_t z = ((real_t)i2 + 0.5) * hgrid_[2];

                data_[m_idx(i0, i1, i2)] = M_PI * pow(x, 4) + M_PI_2 * pow(x, 2) - M_PI_4;
            }
        }
    }

    // get the boundary on the x_left
    PhysBlock*     phys_left = new PhysBlock(0, block_);
    lid_t          fstart[3] = {0, 0, 0};
    EvenBoundary_4 bc;
    bc(phys_left->iface(), fstart, hgrid_, 0.0, phys_left, data_);

    // check the result for LEFT
    for (int i2 = -M_GS; i2 < M_N + M_GS; i2++) {
        for (int i1 = -M_GS; i1 < M_N + M_GS; i1++) {
            for (int i0 = -M_GS; i0 < 0; i0++) {
                real_t x = ((real_t)i0 + 0.5) * hgrid_[0];
                real_t y = ((real_t)i1 + 0.5) * hgrid_[1];
                real_t z = ((real_t)i2 + 0.5) * hgrid_[2];

                EXPECT_NEAR(data_[m_idx(i0, i1, i2)], M_PI * pow(x, 4) + M_PI_2 * pow(x, 2) - M_PI_4, DOUBLE_TOL) << "failed for " << i0 << " " << i1 << " " << i2;
            }
        }
    }

    delete (phys_left);
}

//==============================================================================================================================
TEST_F(valid_Boundary, bc_even_right) {
    // fill the source
    for (int i2 = -M_GS; i2 < M_N + M_GS; i2++) {
        for (int i1 = -M_GS; i1 < M_N + M_GS; i1++) {
            for (int i0 = 0; i0 < M_N; i0++) {
                real_t x = ((real_t)i0 + 0.5) * hgrid_[0] - 1.0;
                real_t y = ((real_t)i1 + 0.5) * hgrid_[1] - 1.0;
                real_t z = ((real_t)i2 + 0.5) * hgrid_[2] - 1.0;

                data_[m_idx(i0, i1, i2)] = M_PI * pow(x, 4) + M_PI_2 * pow(x, 2) + 1.0;
            }
        }
    }

    // get the boundary on the x_right
    PhysBlock*     phys_right = new PhysBlock(1, block_);
    lid_t          fstart[3]  = {M_N - 1, 0, 0};
    EvenBoundary_4 bc;
    // the flux is 4 pi x^3 + pi x -> 5 pi at the interface
    bc(phys_right->iface(), fstart, hgrid_, 0.0, phys_right, data_);

    // check the result for RIGHT
    for (int i2 = -M_GS; i2 < M_N + M_GS; i2++) {
        for (int i1 = -M_GS; i1 < M_N + M_GS; i1++) {
            for (int i0 = M_N; i0 < M_N + M_GS; i0++) {
                real_t x = ((real_t)i0 + 0.5) * hgrid_[0] - 1.0;
                real_t y = ((real_t)i1 + 0.5) * hgrid_[1] - 1.0;
                real_t z = ((real_t)i2 + 0.5) * hgrid_[2] - 1.0;

                EXPECT_NEAR(data_[m_idx(i0, i1, i2)], M_PI * pow(x, 4) + M_PI_2 * pow(x, 2) + 1.0, DOUBLE_TOL) << "failed for " << i0 << " " << i1 << " " << i2;
            }
        }
    }

    delete (phys_right);
}

//==============================================================================================================================
TEST_F(valid_Boundary, bc_odd_left) {
    // fill the source
    for (int i2 = -M_GS; i2 < M_N + M_GS; i2++) {
        for (int i1 = 0; i1 < M_N; i1++) {
            for (int i0 = -M_GS; i0 < M_N + M_GS; i0++) {
                real_t x = ((real_t)i0 + 0.5) * hgrid_[0];
                real_t y = ((real_t)i1 + 0.5) * hgrid_[1];
                real_t z = ((real_t)i2 + 0.5) * hgrid_[2];

                data_[m_idx(i0, i1, i2)] = M_PI * pow(y, 5) + M_PI_2 * pow(y, 3) - M_PI_4 * pow(y, 1);
            }
        }
    }

    // get the boundary on the x_left
    PhysBlock*    phys_left = new PhysBlock(2, block_);
    lid_t         fstart[3] = {0, 0, 0};
    OddBoundary_4 bc;
    bc(phys_left->iface(), fstart, hgrid_, 0.0, phys_left, data_);

    // check the result for LEFT
    for (int i2 = -M_GS; i2 < M_N + M_GS; i2++) {
        for (int i1 = -M_GS; i1 < M_N; i1++) {
            for (int i0 = -M_GS; i0 < M_N + M_GS; i0++) {
                real_t x = ((real_t)i0 + 0.5) * hgrid_[0];
                real_t y = ((real_t)i1 + 0.5) * hgrid_[1];
                real_t z = ((real_t)i2 + 0.5) * hgrid_[2];

                EXPECT_NEAR(data_[m_idx(i0, i1, i2)], M_PI * pow(y, 5) + M_PI_2 * pow(y, 3) - M_PI_4 * pow(y, 1), DOUBLE_TOL) << "failed for " << i0 << " " << i1 << " " << i2;
            }
        }
    }

    delete (phys_left);
}

//==============================================================================================================================
TEST_F(valid_Boundary, bc_odd_right) {
    // fill the source
    for (int i2 = -M_GS; i2 < M_N + M_GS; i2++) {
        for (int i1 = 0; i1 < M_N; i1++) {
            for (int i0 = -M_GS; i0 < M_N + M_GS; i0++) {
                real_t x = ((real_t)i0 + 0.5) * hgrid_[0] - 1.0;
                real_t y = ((real_t)i1 + 0.5) * hgrid_[1] - 1.0;
                real_t z = ((real_t)i2 + 0.5) * hgrid_[2] - 1.0;

                data_[m_idx(i0, i1, i2)] = M_PI * pow(y, 5) + M_PI_2 * pow(y, 3) - M_PI_4 * pow(y, 1);
            }
        }
    }

    // get the boundary on the x_right
    PhysBlock*    phys_right = new PhysBlock(3, block_);
    lid_t         fstart[3]  = {0, M_N - 1, 0};
    OddBoundary_4 bc;
    // the flux is 4 pi x^3 + pi x -> 5 pi at the interface
    bc(phys_right->iface(), fstart, hgrid_, 0.0, phys_right, data_);

    // check the result for RIGHT
    for (int i2 = -M_GS; i2 < M_N + M_GS; i2++) {
        for (int i1 = 0; i1 < M_N + M_GS; i1++) {
            for (int i0 = -M_GS; i0 < M_N + M_GS; i0++) {
                real_t x = ((real_t)i0 + 0.5) * hgrid_[0] - 1.0;
                real_t y = ((real_t)i1 + 0.5) * hgrid_[1] - 1.0;
                real_t z = ((real_t)i2 + 0.5) * hgrid_[2] - 1.0;

                EXPECT_NEAR(data_[m_idx(i0, i1, i2)], M_PI * pow(y, 5) + M_PI_2 * pow(y, 3) - M_PI_4 * pow(y, 1), DOUBLE_TOL) << "failed for " << i0 << " " << i1 << " " << i2;
            }
        }
    }

    delete (phys_right);
}

//==============================================================================================================================
TEST_F(valid_Boundary, bc_extrap_4_left) {
    // fill the source
    for (int i2 = 0; i2 < M_N; i2++) {
        for (int i1 = -M_GS; i1 < M_N + M_GS; i1++) {
            for (int i0 = -M_GS; i0 < M_N + M_GS; i0++) {
                real_t x = ((real_t)i0 + 0.5) * hgrid_[0];
                real_t y = ((real_t)i1 + 0.5) * hgrid_[1];
                real_t z = ((real_t)i2 + 0.5) * hgrid_[2];

                data_[m_idx(i0, i1, i2)] = M_PI * pow(z, 3) + M_PI_2 * pow(z, 2) - M_PI_4 * pow(z, 1) + M_SQRT2;
            }
        }
    }

    // get the boundary on the x_left
    PhysBlock*       phys_left = new PhysBlock(4, block_);
    lid_t            fstart[3] = {0, 0, 0};
    ExtrapBoundary_4 bc;
    bc(phys_left->iface(), fstart, hgrid_, 0.0, phys_left, data_);

    // check the result for LEFT
    for (int i2 = -M_GS; i2 < M_N; i2++) {
        for (int i1 = -M_GS; i1 < M_N + M_GS; i1++) {
            for (int i0 = -M_GS; i0 < M_N + M_GS; i0++) {
                real_t x = ((real_t)i0 + 0.5) * hgrid_[0];
                real_t y = ((real_t)i1 + 0.5) * hgrid_[1];
                real_t z = ((real_t)i2 + 0.5) * hgrid_[2];

                EXPECT_NEAR(data_[m_idx(i0, i1, i2)], M_PI * pow(z, 3) + M_PI_2 * pow(z, 2) - M_PI_4 * pow(z, 1) + M_SQRT2, DOUBLE_TOL) << "failed for " << i0 << " " << i1 << " " << i2;
            }
        }
    }

    delete (phys_left);
}
TEST_F(valid_Boundary, bc_extrap_4_right) {
    // fill the source
    for (int i2 = 0; i2 < M_N; i2++) {
        for (int i1 = -M_GS; i1 < M_N + M_GS; i1++) {
            for (int i0 = -M_GS; i0 < M_N + M_GS; i0++) {
                real_t x = ((real_t)i0 + 0.5) * hgrid_[0] - 1.0;
                real_t y = ((real_t)i1 + 0.5) * hgrid_[1] - 1.0;
                real_t z = ((real_t)i2 + 0.5) * hgrid_[2] - 1.0;

                data_[m_idx(i0, i1, i2)] = M_PI * pow(z, 3) + M_PI_2 * pow(z, 2) - M_PI_4 * pow(z, 1) + M_SQRT2;
            }
        }
    }

    // get the boundary on the x_right
    PhysBlock*       phys_right = new PhysBlock(5, block_);
    lid_t            fstart[3]  = {0, 0, M_N - 1};
    ExtrapBoundary_4 bc;
    // the flux is 4 pi x^3 + pi x -> 5 pi at the interface
    bc(phys_right->iface(), fstart, hgrid_, 0.0, phys_right, data_);

    // check the result for RIGHT
    for (int i2 = 0; i2 < M_N + M_GS; i2++) {
        for (int i1 = -M_GS; i1 < M_N + M_GS; i1++) {
            for (int i0 = -M_GS; i0 < M_N + M_GS; i0++) {
                real_t x = ((real_t)i0 + 0.5) * hgrid_[0] - 1.0;
                real_t y = ((real_t)i1 + 0.5) * hgrid_[1] - 1.0;
                real_t z = ((real_t)i2 + 0.5) * hgrid_[2] - 1.0;

                EXPECT_NEAR(data_[m_idx(i0, i1, i2)], M_PI * pow(z, 3) + M_PI_2 * pow(z, 2) - M_PI_4 * pow(z, 1) + M_SQRT2, DOUBLE_TOL) << "failed for " << i0 << " " << i1 << " " << i2;
            }
        }
    }

    delete (phys_right);
}
//==============================================================================================================================
TEST_F(valid_Boundary, bc_extrap_3_left) {
    // fill the source
    for (int i2 = 0; i2 < M_N; i2++) {
        for (int i1 = -M_GS; i1 < M_N + M_GS; i1++) {
            for (int i0 = -M_GS; i0 < M_N + M_GS; i0++) {
                real_t x = ((real_t)i0 + 0.5) * hgrid_[0];
                real_t y = ((real_t)i1 + 0.5) * hgrid_[1];
                real_t z = ((real_t)i2 + 0.5) * hgrid_[2];

                data_[m_idx(i0, i1, i2)] = M_PI_2 * pow(z, 2) - M_PI_4 * pow(z, 1) + M_SQRT2;
            }
        }
    }

    // get the boundary on the x_left
    PhysBlock*       phys_left = new PhysBlock(4, block_);
    lid_t            fstart[3] = {0, 0, 0};
    ExtrapBoundary_3 bc;
    bc(phys_left->iface(), fstart, hgrid_, 0.0, phys_left, data_);

    // check the result for LEFT
    for (int i2 = -M_GS; i2 < M_N; i2++) {
        for (int i1 = -M_GS; i1 < M_N + M_GS; i1++) {
            for (int i0 = -M_GS; i0 < M_N + M_GS; i0++) {
                real_t x = ((real_t)i0 + 0.5) * hgrid_[0];
                real_t y = ((real_t)i1 + 0.5) * hgrid_[1];
                real_t z = ((real_t)i2 + 0.5) * hgrid_[2];

                EXPECT_NEAR(data_[m_idx(i0, i1, i2)], M_PI_2 * pow(z, 2) - M_PI_4 * pow(z, 1) + M_SQRT2, DOUBLE_TOL) << "failed for " << i0 << " " << i1 << " " << i2;
            }
        }
    }

    delete (phys_left);
}
TEST_F(valid_Boundary, bc_extrap_3_right) {
    // fill the source
    for (int i2 = 0; i2 < M_N; i2++) {
        for (int i1 = -M_GS; i1 < M_N + M_GS; i1++) {
            for (int i0 = -M_GS; i0 < M_N + M_GS; i0++) {
                real_t x = ((real_t)i0 + 0.5) * hgrid_[0] - 1.0;
                real_t y = ((real_t)i1 + 0.5) * hgrid_[1] - 1.0;
                real_t z = ((real_t)i2 + 0.5) * hgrid_[2] - 1.0;

                data_[m_idx(i0, i1, i2)] = M_PI_2 * pow(z, 2) - M_PI_4 * pow(z, 1) + M_SQRT2;
            }
        }
    }

    // get the boundary on the x_right
    PhysBlock*       phys_right = new PhysBlock(5, block_);
    lid_t            fstart[3]  = {0, 0, M_N - 1};
    ExtrapBoundary_3 bc;
    // the flux is 4 pi x^3 + pi x -> 5 pi at the interface
    bc(phys_right->iface(), fstart, hgrid_, 0.0, phys_right, data_);

    // check the result for RIGHT
    for (int i2 = 0; i2 < M_N + M_GS; i2++) {
        for (int i1 = -M_GS; i1 < M_N + M_GS; i1++) {
            for (int i0 = -M_GS; i0 < M_N + M_GS; i0++) {
                real_t x = ((real_t)i0 + 0.5) * hgrid_[0] - 1.0;
                real_t y = ((real_t)i1 + 0.5) * hgrid_[1] - 1.0;
                real_t z = ((real_t)i2 + 0.5) * hgrid_[2] - 1.0;

                EXPECT_NEAR(data_[m_idx(i0, i1, i2)], M_PI_2 * pow(z, 2) - M_PI_4 * pow(z, 1) + M_SQRT2, DOUBLE_TOL) << "failed for " << i0 << " " << i1 << " " << i2;
            }
        }
    }

    delete (phys_right);
}

//==============================================================================================================================
TEST_F(valid_Boundary, bc_extrap_5_left) {
    // fill the source
    for (int i2 = 0; i2 < M_N; i2++) {
        for (int i1 = -M_GS; i1 < M_N + M_GS; i1++) {
            for (int i0 = -M_GS; i0 < M_N + M_GS; i0++) {
                real_t x = ((real_t)i0 + 0.5) * hgrid_[0];
                real_t y = ((real_t)i1 + 0.5) * hgrid_[1];
                real_t z = ((real_t)i2 + 0.5) * hgrid_[2];

                data_[m_idx(i0, i1, i2)] = M_PI * pow(z, 4) + M_PI_2 * pow(z, 3) - M_PI_4 * pow(z, 2) + M_SQRT2;
            }
        }
    }

    // get the boundary on the x_left
    PhysBlock*       phys_left = new PhysBlock(4, block_);
    lid_t            fstart[3] = {0, 0, 0};
    ExtrapBoundary_5 bc;
    bc(phys_left->iface(), fstart, hgrid_, 0.0, phys_left, data_);

    // check the result for LEFT
    for (int i2 = -M_GS; i2 < M_N; i2++) {
        for (int i1 = -M_GS; i1 < M_N + M_GS; i1++) {
            for (int i0 = -M_GS; i0 < M_N + M_GS; i0++) {
                real_t x = ((real_t)i0 + 0.5) * hgrid_[0];
                real_t y = ((real_t)i1 + 0.5) * hgrid_[1];
                real_t z = ((real_t)i2 + 0.5) * hgrid_[2];

                ASSERT_NEAR(data_[m_idx(i0, i1, i2)], M_PI * pow(z, 4) + M_PI_2 * pow(z, 3) - M_PI_4 * pow(z, 2) + M_SQRT2, DOUBLE_TOL) << "failed for " << i0 << " " << i1 << " " << i2;
            }
        }
    }

    delete (phys_left);
}
TEST_F(valid_Boundary, bc_extrap_5_right) {
    // fill the source
    for (int i2 = 0; i2 < M_N; i2++) {
        for (int i1 = -M_GS; i1 < M_N + M_GS; i1++) {
            for (int i0 = -M_GS; i0 < M_N + M_GS; i0++) {
                real_t x = ((real_t)i0 + 0.5) * hgrid_[0] - 1.0;
                real_t y = ((real_t)i1 + 0.5) * hgrid_[1] - 1.0;
                real_t z = ((real_t)i2 + 0.5) * hgrid_[2] - 1.0;

                data_[m_idx(i0, i1, i2)] = M_PI * pow(z, 4) + M_PI_2 * pow(z, 3) - M_PI_4 * pow(z, 2) + M_SQRT2;
            }
        }
    }

    // get the boundary on the x_right
    PhysBlock*       phys_right = new PhysBlock(5, block_);
    lid_t            fstart[3]  = {0, 0, M_N - 1};
    ExtrapBoundary_5 bc;
    // the flux is 4 pi x^3 + pi x -> 5 pi at the interface
    bc(phys_right->iface(), fstart, hgrid_, 0.0, phys_right, data_);

    // check the result for RIGHT
    for (int i2 = 0; i2 < M_N + M_GS; i2++) {
        for (int i1 = -M_GS; i1 < M_N + M_GS; i1++) {
            for (int i0 = -M_GS; i0 < M_N + M_GS; i0++) {
                real_t x = ((real_t)i0 + 0.5) * hgrid_[0] - 1.0;
                real_t y = ((real_t)i1 + 0.5) * hgrid_[1] - 1.0;
                real_t z = ((real_t)i2 + 0.5) * hgrid_[2] - 1.0;

                ASSERT_NEAR(data_[m_idx(i0, i1, i2)],M_PI * pow(z, 4) + M_PI_2 * pow(z, 3) - M_PI_4 * pow(z, 2) + M_SQRT2, DOUBLE_TOL) << "failed for " << i0 << " " << i1 << " " << i2;
            }
        }
    }

    delete (phys_right);
}