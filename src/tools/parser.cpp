#include "parser.hpp"

#include <argp.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include <string>

#include "core/macros.hpp"
#include "core/types.hpp"
#include "navier_stokes.hpp"

using std::string;

static int count_list(const char* arg) {
    int         count = 0;
    const char* tmp   = strchr(arg, ',');
    while (tmp != NULL) {
        count++;
        tmp = strchr(tmp + 1, ',');
    }
    return count + 1;
}

static error_t atoc_list(const int length, char* arg, char* list) {
    const int count = count_list(arg);
    if (count != length) {
        return ARGP_ERR_UNKNOWN;
    } else {
        char* num = strtok(arg, ",");
        for (int id = 0; id < length; id++) {
            int temp = atoi(num);
            m_assert(std::numeric_limits<char>::max() > temp, "the number read is too big to be cast in a char");
            list[id] = (char)(temp);
            num      = strtok(NULL, ",");
        }
    }
    return 0;
} 

static error_t atoi_list(const int length, char* arg, int* list) {
    const int count = count_list(arg);
    if (count != length) {
        return ARGP_ERR_UNKNOWN;
    } else {
        char* num = strtok(arg, ",");
        for (int id = 0; id < length; id++) {
            list[id] = atoi(num);
            num      = strtok(NULL, ",");
        }
    }
    return 0;
}

static error_t atob_list(const int length, char* arg, bool* list) {
    const int count = count_list(arg);
    if (count != length) {
        return ARGP_ERR_UNKNOWN;
    } else {
        char* num = strtok(arg, ",");
        for (int id = 0; id < length; id++) {
            int mybool = atoi(num);
            if (!(mybool == 0 || mybool == 1)) {
                return ARGP_ERR_UNKNOWN;
            }
            list[id] = (bool)mybool;
            num      = strtok(NULL, ",");
        }
    }
    return 0;
}

static error_t atof_list(const int length, char* arg, real_t* list) {
    const int count = count_list(arg);
    if (count != length) {
        return ARGP_ERR_UNKNOWN;
    } else {
        char* num = strtok(arg, ",");
        for (int id = 0; id < length; id++) {
            list[id] = atof(num);
            num      = strtok(NULL, ",");
        }
    }
    return 0;
}

static char doc[] = "MURPHY - a MUltiResolution multiPHYsics framework.";

static struct argp_option options[] = {
    /* general parameters */
    {0, 0, 0, OPTION_DOC, "Options:", 1},
    {"profile", 'p', 0, OPTION_ARG_OPTIONAL, "enables the program profiling", 1},

    /* general param */
    {0, 0, 0, OPTION_DOC, "Grid-related parameters:", 2},
    {"dom", 2001, "d_x,d_y,d_z", 0, "gives the dimension of the domain (integers: d_x,d_y,d_z)"},
    {"ilevel", 2002, "level", 0, "the initialization level (integer: num)"},
    {"periodic", 2003, "p_x,p_y,p_z", 0, "periodicity of the domain"},
    {"patch", 2004, "o_x,o_y,o_z,l_x,l_y,l_z,lvl", 0, "patch from (o_x,o_y,o_z) and of length (l_x,l_y,l_z) at level (lvl)"},
    {"rtol", 2005, "tol", 0, "refinement tolerance"},
    {"ctol", 2006, "tol", 0, "coarsening tolerance"},
    {"no-adapt", 2007, 0, OPTION_ARG_OPTIONAL, "disable the mesh adaptation"},
    {"level-min", 2008, "level", 0, "the minimum level on the grid (integer: num)"},
    {"level-max", 2009, "level", 0, "the maximum level on the grid (integer: num)"},

    /* general parameters */
    {0, 0, 0, OPTION_DOC, "Other parameters:", 3},
    {"reynolds", 3001, "double", 0, "the Reynolds number"},
    {"vr-normal", 3002, "dir", 0, "vortex ring normal"},
    {"vr-center", 3003, "c_x,c_y,c_z", 0, "vortex ring center"},
    {"vr-radius", 3004, "rad", 0, "vortex ring radius"},
    {"vr-sigma", 3005, "sigma", 0, "vortex ring sigma"},
    {"compute-error", 3006, "bool", 0, "compute the error when running diagnostics"},
    {"dump-error", 3007, "bool", 0, "dump the error when running diagnostics and if computer error is true"},
    {"dump-detail", 3008, "bool", 0, "dump the detail coefficients when running diagnostics"},
    {"iter-max", 3009, "int", 0, "maximum of RK3 iterations"},
    {"iter-diag", 3010, "int", 0, "run the diagnostics every x iterations"},
    {"iter-adapt", 3011, "int", 0, "adapt the grid every x iterations"},
    {"iter-dump", 3012, "int", 0, "dump the field every x iterations"},
    {"no-weno", 3013, 0, OPTION_ARG_OPTIONAL, "disable the weno adaptation and uses a regular upwind-downwind stencil"},
    {"grid-on-sol", 3014, 0, OPTION_ARG_OPTIONAL, "adapt the grid based on the solution"},
    {"weno-5", 3015, 0, OPTION_ARG_OPTIONAL, "uses the 5th order WENO stencil"},
    {"eps-start", 3016, "value", 0, "start epsilon"},
    {"delta-eps", 3017, "value", 0, "factor from one epsilon to another"},

    /* client choice parameters */
    {0, 0, 0, OPTION_DOC, "Available clients:", 4},
    // navier-stokes
    {"navier-stokes", 4000, 0, OPTION_ARG_OPTIONAL, "Navier-Stokes testcase"},
    {"ns", 0, 0, OPTION_ALIAS, 0},
    // abc flow
    {"abc", 4001, 0, OPTION_ARG_OPTIONAL, "ABC-flow testcase"},
    // simple advection
    {"simple-advection", 4002, 0, OPTION_ARG_OPTIONAL, "simple advection testcase"},
    {"sadv", 0, 0, OPTION_ALIAS, 0},
    // simple advection
    {"epsilon-test", 4003, 0, OPTION_ARG_OPTIONAL, "Epsilon test"},
    {"eps", 0, 0, OPTION_ALIAS, 0},

    /* help */
    {0, 0, 0, OPTION_DOC, "Help:", -1},
    {0}
    };

static error_t parse_opt(int key, char* arg, struct argp_state* state) {
    ParserArguments* arguments = reinterpret_cast<ParserArguments*>(state->input);

    switch (key) {
        //................................................
        case 'p': { /* profiler */
            arguments->profile = true;
            m_log("profiler set to true");
            return 0;
        }
        //................................................
        case 2001: { /* domain aspect ratio */
            int*    length = arguments->length;
            error_t err    = atoi_list(3, arg, length);
            m_log("domain length: %d %d %d", length[0], length[1], length[2]);
            return err;
        }
        case 2002: { /* level */
            char*    lvl = &arguments->init_lvl;
            error_t err = atoc_list(1, arg, lvl);
            m_log("init level: %d", lvl[0]);
            return err;
        }
        case 2003: { /* periodic */
            bool*   period = arguments->period;
            error_t err    = atob_list(3, arg, period);
            m_log("periodicity: %d %d %d", period[0], period[1], period[2]);
            return err;
        }
        case 2004: { /* patches */
            real_t  data[7];
            error_t err   = atof_list(7, arg, data);
            lid_t   level = (lid_t)(data[6]);
            arguments->patch.push_back(Patch(data, data + 3, level));
            m_log("patch: level %d starting (%f %f %f) of length (%f %f %f)", level, data[0], data[1], data[2], data[3], data[4], data[5]);
            return err;
        }
        case 2005: { /* rtol */
            real_t* tol = &arguments->refine_tol;
            error_t err = atof_list(1, arg, tol);
            m_log("refinement tolerance: %f", tol[0]);
            return err;
        }
        case 2006: { /* ctol */
            real_t* tol = &arguments->coarsen_tol;
            error_t err = atof_list(1, arg, tol);
            m_log("coarsening tolerance: %f", tol[0]);
            return err;
        }
        case 2007: { /* no adapt */
            arguments->no_adapt = true;
            m_log("no_adaptation");
            return 0;
        }
        case 2008: { /* level_min */
            level_t* lvl = &arguments->level_min;
            error_t  err = atoc_list(1, arg, lvl);
            m_log("level min: %d", lvl[0]);
            return err;
        }
        case 2009: { /* level_max */
            level_t* lvl = &arguments->level_max;
            error_t  err = atoc_list(1, arg, lvl);
            m_log("level max: %d", lvl[0]);
            return err;
        }
        //................................................
        case 3001: { /* Reynolds */
            double* reynolds = &arguments->reynolds;
            error_t err      = atof_list(1, arg, reynolds);
            m_log("Reynolds number: %f", reynolds[0]);
            return err;
        }
        case 3002: { /* VR normal */
            int* normal = &arguments->vr_normal;
            error_t err      = atoi_list(1, arg, normal);
            m_log("vortex-ring normal: %d", normal[0]);
            return err;
        }
        case 3003: { /* VR center */
            real_t* center = arguments->vr_center;
            error_t err      = atof_list(3, arg, center);
            m_log("vortex-ring center: %f %f %f", center[0], center[1], center[2]);
            return err;
        }
        case 3004: { /* VR radius */
            real_t* rad = &arguments->vr_radius;
            error_t err      = atof_list(1, arg, rad);
            m_log("vortex-ring radius: %f", rad[0]);
            return err;
        }
        case 3005: { /* VR sigma */
            real_t* sigma = &arguments->vr_sigma;
            error_t err      = atof_list(1, arg, sigma);
            m_log("vortex-ring sigma: %f", sigma[0]);
            return err;
        }
        case 3006: { /* compute error*/
            bool*   mybool = &arguments->compute_error;
            error_t err    = atob_list(1, arg, mybool);
            m_log("compute error: %d", mybool[0]);
            return err;
        }
        case 3007: { /* dump error*/
            bool*   mybool = &arguments->dump_error;
            error_t err    = atob_list(1, arg, mybool);
            m_log("dump error: %d", mybool[0]);
            return err;
        }
        case 3008: { /* dump details*/
            bool*   mybool = &arguments->dump_detail;
            error_t err    = atob_list(1, arg, mybool);
            m_log("dump detail: %d", mybool[0]);
            return err;
        }
        case 3009: { /* iter max*/
            int*    myint = &arguments->iter_max;
            error_t err   = atoi_list(1, arg, myint);
            m_log("iter max: %d", myint[0]);
            return err;
        }
        case 3010: { /* iter diag*/
            int*    myint = &arguments->iter_diag;
            error_t err   = atoi_list(1, arg, myint);
            m_log("iter diag: %d", myint[0]);
            return err;
        }
        case 3011: { /* iter adapt*/
            int*    myint = &arguments->iter_adapt;
            error_t err   = atoi_list(1, arg, myint);
            m_log("iter adapt: %d", myint[0]);
            return err;
        }
        case 3012: { /* iter dump*/
            int*    myint = &arguments->iter_dump;
            error_t err   = atoi_list(1, arg, myint);
            m_log("iter dump: %d", myint[0]);
            return err;
        }
        case 3013: { /* no-weno */
            arguments->no_weno = true;
            m_log("no-weno");
            return 0;
        }
        case 3014: { /* grid-on-sol */
            arguments->grid_on_sol = true;
            m_log("grid based on the solution");
            return 0;
        }
        case 3015: { /* weno_5 */
            arguments->weno_5 = true;
            m_log("weno");
            return 0;
        }
        case 3016: { /* rtol */
            real_t* tol = &arguments->eps_start;
            error_t err = atof_list(1, arg, tol);
            m_log("epsilon start: %f", tol[0]);
            return err;
        }
        case 3017: { /* rtol */
            real_t* tol = &arguments->delta_eps;
            error_t err = atof_list(1, arg, tol);
            m_log("delta epsilon: %f", tol[0]);
            return err;
        }
        //................................................
        case 4000: { /* Navier-Stockes */
            m_log("Navier-Stokes testcase selected");
            arguments->do_navier_stokes = true;
            return 0;
        }
        case 4001: { /* ABC flow */
            m_log("ABC-flow testcase selected");
            arguments->do_abc_flow = true;
            return 0;
        }
        case 4002: { /* simple advection */
            m_log("Simple-advection testcase selected");
            arguments->do_simple_adv = true;
            return 0;
        }
        case 4003: { /* simple advection */
            m_log("Epsilon testcase selected");
            arguments->do_epsilon_test = true;
            return 0;
        }
        
        default:
            return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

struct argp argp = {options, parse_opt, 0, doc};

void ParseArgument(int argc, char** argv, ParserArguments* arguments) {
    argp_parse(&argp, argc, argv, ARGP_NO_EXIT, 0, arguments);
}