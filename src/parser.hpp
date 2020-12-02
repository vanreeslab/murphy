#ifndef SRC_PARSER_HPP_
#define SRC_PARSER_HPP_

#include <argp.h>

#include <list>

#include "defs.hpp"
#include "patch.hpp"

class ParserArguments {
   public:
    int              init_lvl  = 0;                      //!< the level at which the grid is initialized
    bool             period[3] = {false, false, false};  //!< the periodicity of the domain
    lid_t            length[3] = {1, 1, 1};              //!< the aspect ratio of the domain
    std::list<Patch> patch;                              //!< list of imposed level regions for the initialization

    bool   profile          = false;
    bool   do_navier_stokes = false;
    real_t reynolds         = 0.0;

    real_t refine_tol  = 1e-1;
    real_t coarsen_tol = 1e-4;

    int    vr_normal    = 2;
    real_t vr_radius    = 0.25;
    real_t vr_sigma     = 0.025;
    real_t vr_center[3] = {0.5, 0.5, 0.5};
};

void ParseArgument(int argc, char** argv, ParserArguments* arguments);

#endif
