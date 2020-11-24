#include "parser.hpp"

#include <argp.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include <string>

#include "defs.hpp"
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
    {"profile", 1001, 0, OPTION_ARG_OPTIONAL, "enables the program profiling", 1},

    /* general param */
    {0, 0, 0, OPTION_DOC, "General parameters:", 2},
    {"dom", 2001, "d_x,d_y,d_z", 0, "gives the dimension of the domain (integers: d_x,d_y,d_z)"},
    {"ilevel", 2002, "level", 0, "the initialization level (integer: num)"},
    {"periodic", 2003, "p_x,p_y,p_z", 0, "periodicity of the domain"},
    {"patch", 2004, "o_x,o_y,o_z,l_x,l_y,l_z,lvl", 0, "patch from (o_x,o_y,o_z) and of length (l_x,l_y,l_z) at level (lvl)"},
    {"rtol", 2005, "tol", 0, "refinement tolerance"},
    {"ctol", 2006, "tol", 0, "coarsening tolerance"},

    /* client choice parameters */
    {0, 0, 0, OPTION_DOC, "Available clients:", 3},
    {"navier-stokes", 3000, 0, OPTION_ARG_OPTIONAL, "Navier-Stokes testcase"},
    {"ns", 0, 0, OPTION_ALIAS, 0},

    /* Navier-Stokes */
    {0, 0, 0, OPTION_DOC, "Navier-Stokes parameters:", 4},
    {"reynolds", 3001, "NUM", 0, "the Reynolds number"},
    {"vr-normal", 3002, "dir", 0, "vortex ring normal"},
    {"vr-center", 3003, "c_x,c_y,c_z", 0, "vortex ring center"},
    {"vr-radius", 3004, "rad", 0, "vortex ring radius"},
    {"vr-sigma", 3005, "sigma", 0, "vortex ring sigma"},

    /* help */
    {0, 0, 0, OPTION_DOC, "Help:", -1},
    {0}};

static error_t parse_opt(int key, char* arg, struct argp_state* state) {
    ParserArguments* arguments = reinterpret_cast<ParserArguments*>(state->input);

    switch (key) {
        //................................................
        case 1001: { /* profiler */
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
            int*    lvl = &arguments->init_lvl;
            error_t err = atoi_list(1, arg, lvl);
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
            m_log("refinement tolerance: %d", tol[0]);
            return err;
        }
        case 2006: { /* ctol */
            real_t* tol = &arguments->coarsen_tol;
            error_t err = atof_list(1, arg, tol);
            m_log("coarsening tolerance: %d", tol[0]);
            return err;
        }
        //................................................
        case 3000: { /* Navier-Stockes */
            m_log("Navier-Stockes");
            arguments->do_navier_stokes = true;
            return 0;
        }
        case 3001: { /* Reynolds */
            double* reynolds = &arguments->reynolds;
            error_t err      = atof_list(1, arg, reynolds);
            m_log("Reynolds number: %f", reynolds[0]);
            return err;
        }
        default:
            return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

static struct argp argp = {options, parse_opt, 0, doc};

void ParseArgument(int argc, char** argv, ParserArguments* arguments) {
    argp_parse(&argp, argc, argv, 0, 0, arguments);
}