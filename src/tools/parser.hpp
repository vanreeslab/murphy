#ifndef SRC_PARSER_HPP_
#define SRC_PARSER_HPP_

#include <argp.h>

#include <list>

#include "core/macros.hpp"
#include "core/types.hpp"
#include "patch.hpp"

struct ParserArguments {
    level_t          init_lvl  = 0;                      //!< the level at which the grid is initialized
    bool             period[3] = {false, false, false};  //!< the periodicity of the domain
    lid_t            length[3] = {1, 1, 1};              //!< the aspect ratio of the domain
    std::list<Patch> patch;                              //!< list of imposed level regions for the initialization

    bool profile          = false;
    bool do_navier_stokes = false;
    bool do_abc_flow      = false;
    bool do_simple_adv    = false;
    bool do_epsilon_test  = false;
    bool do_debug_lifting = false;
    bool do_conv_weno     = false;
    bool do_2lvl_weno     = false;
    bool do_weak_scal     = false;

    real_t reynolds    = 0.0;
    real_t refine_tol  = 1e-2;
    real_t coarsen_tol = 1e-4;
    bool   optimal_tol = false;

    real_t time_start = 0.0;
    real_t time_final = 0.5;

    real_t eps_start = 1.0;
    real_t delta_eps = 0.1;

    real_t cfl_max = 100.0;

    level_t level_min = 0;
    level_t level_max = P8EST_QMAXLEVEL;

    bool no_adapt = false;

    int    vr_normal    = 2;
    real_t vr_radius    = 0.25;
    real_t vr_sigma     = 0.025;
    real_t vr_center[3] = {0.5, 0.5, 0.5};

    bool   dump_error    = false;
    bool   dump_detail   = false;
    bool   compute_error = true;
    iter_t iter_max      = 100;
    iter_t iter_diag     = 1;
    iter_t iter_adapt    = 1;
    iter_t iter_dump     = 1;
    // bool   no_weno       = false;
    // bool   weno_5        = false;
    bool fix_weno    = false;
    int  weno        = 3;
    bool grid_on_sol = false;
};

void ParseArgument(int argc, char** argv, ParserArguments* arguments);

#endif
