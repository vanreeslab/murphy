#include <argp.h>

#include "navier_stokes.hpp"
#include "parser.hpp"

static struct argp_option options[] = {
    {"domain", 'd', "d_x,d_y,d_z", 0, "dimension of the domain in x,y,z (int)"},
    {"reynolds", 'r', "num", 0, "the Reynolds number (double)"},
    {0}};

typedef struct NavierStokesParam {
    real_t Re_   = 0.0;        //!< the Reynolds number
    lid_t  L_[3] = {1, 1, 5};  //!< the dimension of the domain
} NavierStokesParam;

static error_t parse_opt(int key, char* arg, struct argp_state* state) {
    NavierStokesParam* arguments = reinterpret_cast<NavierStokesParam*>(state->input);

    switch (key) {
        case 'r': {
            real_t* Re  = &arguments->Re_;
            error_t err = atof_list(1, arg, Re);
            m_log("Reynolds: %f", *Re);
            return err;
        }
        case 'd': {
            int*    length = arguments->L_;
            error_t err    = atoi_list(3, arg, length);
            m_log("domain length: %d %d %d", length[0], length[1], length[2]);
            return err;
        }
        default:
            return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

// need to define a global void argp
struct argp extern_ns_argp = {options, 0, 0, 0, 0};          //!< structure used to display the help
struct argp ns_argp        = {options, parse_opt, 0, 0, 0};  //!< structure used to parse the meaningful information

NavierStokes::NavierStokes() {
    //-------------------------------------------------------------------------
    // parse the command line to get info

    //-------------------------------------------------------------------------
}