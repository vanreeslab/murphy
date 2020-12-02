#ifndef SRC_NAVIER_STOKES_HPP_
#define SRC_NAVIER_STOKES_HPP_

#include <string>

#include "defs.hpp"
#include "grid.hpp"
#include "parser.hpp"
#include "prof.hpp"
#include "testcase.hpp"

class NavierStokes : public TestCase {
    real_t nu_ = 0.0;
    real_t u_stream_[3];

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