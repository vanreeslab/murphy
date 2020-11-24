#include "parser.hpp"

#include <argp.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include <string>

#include "murphy.hpp"

// using std::string;

// static char doc[] = "MUPHY - a multiresolution multiphysics framework.";

// static struct argp_option options[] = {
//     // {"periodic", 10001, "p_x,p_y,p_z", 0, "gives the periodicity (boolean: p_x,p_y,p_z)"},
//     // {"domain", 'd', "d_x,d_y,d_z", 0, "gives the dimension of the domain (integers: d_x,d_y,d_z)"},
//     // {"repeat", 'r', "num", 0, "repeat the test num times (integer: num)"},
//     // {"ilevel", 'l', "num", 0, "the initialization level (integer: num)"},
//     // {"patch", 10002, "o_x,o_y,o_z,l_x,l_y,l_z,lvl", 0, "indicate a patch of origin (floats: o_x,o_y,o_z), of length (floats: l_x,l_y,l_z) at level (integer: lvl)"},
//     {0}};

// static error_t parse_opt(int key, char* arg, struct argp_state* state) {
//     parse_arg_t* arguments = reinterpret_cast<parse_arg_t*>(state->input);

//     switch (key) {
//         // case 'r': {
//         //     int*    repeat = &arguments->n_repeat_;
//         //     error_t err    = atoi_list(1, arg, repeat);
//         //     m_log("repeat: %d", repeat[0]);
//         //     return err;
//         // }
//         // case 'l': {
//         //     int*    lvl = &arguments->init_lvl_;
//         //     error_t err = atoi_list(1, arg, lvl);
//         //     m_log("init level: %d", lvl[0]);
//         //     return err;
//         // }
//         // case 'd': {
//         //     int*    length = arguments->length_;
//         //     error_t err    = atoi_list(3, arg, length);
//         //     m_log("domain length: %d %d %d", length[0], length[1], length[2]);
//         //     return err;
//         // }
//         // case 10001: {
//         //     bool*   period = arguments->period_;
//         //     error_t err    = atob_list(3, arg, period);
//         //     m_log("periodicity: %d %d %d", period[0], period[1], period[2]);
//         //     return err;
//         // }
//         // case 10002: {
//         //     real_t  data[7];
//         //     error_t err   = atof_list(7, arg, data);
//         //     lid_t   level = (lid_t)(data[6]);
//         //     arguments->patch_.push_back(Patch(data, data + 3, level));
//         //     m_log("patch: level %d starting (%f %f %f) of length (%f %f %f)", level, data[0], data[1], data[2], data[3], data[4], data[5]);
//         //     return err;
//         // }
//         default:
//             return ARGP_ERR_UNKNOWN;
//     }
//     return 0;
// }

// static struct argp argp = {options, parse_opt, 0, doc};

// void ParseArgument(int argc, char** argv, parse_arg_t* arguments) {
//     argp_parse(&argp, argc, argv, 0, 0, arguments);
// }