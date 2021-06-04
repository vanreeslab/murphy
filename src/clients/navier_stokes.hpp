#ifndef SRC_NAVIER_STOKES_HPP_
#define SRC_NAVIER_STOKES_HPP_

#include <string>

#include "core/macros.hpp"
#include "core/types.hpp"
#include "grid.hpp"
#include "parser.hpp"
#include "prof.hpp"
#include "testcase.hpp"

class NavierStokes : public TestCase {
    bool compute_error_ = false;
    bool dump_error_    = false;

    real_t nu_          = 0.0;
    real_t u_stream_[3] = {0.0, 0.0, 0.0};

    Grid*  grid_ = nullptr;
    Prof*  prof_ = nullptr;
    Field* vort_ = nullptr;

    std::string folder_diag_ = "data";

   public:
    ~NavierStokes();

    void InitParam(ParserArguments* param) override;
    void Run() override;

   protected:
    void Diagnostics(const real_t time, const real_t dt, const lid_t iter);
};

#endif  // SRC_NAVIER_STOKES_HPP_