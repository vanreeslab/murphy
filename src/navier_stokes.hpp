#ifndef SRC_NAVIER_STOKES_HPP_
#define SRC_NAVIER_STOKES_HPP_

#include "defs.hpp"
#include "grid.hpp"
#include "parser.hpp"
#include "prof.hpp"
#include "testcase.hpp"

// static char        ns_doc[] = "Navier-Stokes";
// extern struct argp extern_ns_argp; //!< promise the declaration of a struct argp somewhere

class NavierStokes : public TestCase {

    real_t reynolds_ = 0.0;
    real_t u_stream_[3];

    Grid*  grid_ = nullptr;
    Prof*  prof_ = nullptr;
    Field* vort_ = nullptr;

   public:
    explicit NavierStokes(){};
    ~NavierStokes();

    void InitParam(ParserArguments* param);
    void Run();
};

#endif  // SRC_NAVIER_STOKES_HPP_