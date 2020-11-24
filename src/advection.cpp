// #include "advection.hpp"

// #include <argp.h>

// #include "parser.hpp"

// static struct argp_option options[] = {
//     {"domain", 'd', "d_x,d_y,d_z", 0, "dimension of the domain in x,y,z (int)"},
//     {"periodic", 'p', "p_x,p_y,p_z", 0, "periodicity of the domain in x,y,z (bool)"},
//     {0}};

// typedef struct AdvectionParam {
//     lid_t L_[3]      = {1, 1, 5};  //!< the dimension of the domain
//     bool  period[3] = {0, 0, 0};
// } AdvectionParam;

// static error_t parse_opt(int key, char* arg, struct argp_state* state) {
//     AdvectionParam* arguments = reinterpret_cast<AdvectionParam*>(state->input);

//     switch (key) {
//         case 'd': {
//             int*    length = arguments->L_;
//             error_t err    = atoi_list(3, arg, length);
//             m_log("advection: domain length: %d %d %d", length[0], length[1], length[2]);
//             return err;
//         }
//         case 'p': {
//             bool*   period = arguments->period;
//             error_t err    = atob_list(3, arg, period);
//             m_log("advection: periodicity: %d %d %d", period[0], period[1], period[2]);
//             return err;
//         }
//         default:
//             return ARGP_ERR_UNKNOWN;
//     }
//     return 0;
// }

// struct argp adv_argp        = {options, parse_opt, 0, 0}; //< to parse the option
// struct argp extern_adv_argp = {options, 0, 0, 0}; //!< for the display of the help only

// Advection::Advection() {
//     //-------------------------------------------------------------------------
//     // parse the command line to get info
//     AdvectionParam param;

//     //-------------------------------------------------------------------------
// }