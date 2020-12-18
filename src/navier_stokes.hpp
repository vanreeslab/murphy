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
    bool dump_details_  = false;
    bool compute_error_ = false;
    bool dump_error_    = false;

    real_t nu_          = 0.0;
    real_t u_stream_[3] = {0.0, 0.0, 0.0};
    iter_t iter_max_    = 0;
    iter_t iter_diag_   = 1;
    iter_t iter_adapt_  = 1;

    Grid*  grid_ = nullptr;
    Prof*  prof_ = nullptr;
    Field* vort_ = nullptr;

    std::string folder_diag_ = "data";

   public:
    explicit NavierStokes(){};
    ~NavierStokes();

    void InitParam(ParserArguments* param);
    void Run();
    void Diagnostics(const real_t time, const real_t dt, const lid_t iter);
};

#endif  // SRC_NAVIER_STOKES_HPP_