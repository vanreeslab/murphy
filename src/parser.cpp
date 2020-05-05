#include "parser.hpp"

#include <argp.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include <string>

#include "murphy.hpp"

using std::string;

static char doc[] = "MUPHY - a multiresolution multiphysics framework.";

static struct argp_option options[] = {
    {"periodic", 10001, "p_x,p_y,p_z", 0, "gives the periodicity (boolean: p_x,p_y,p_z)"},
    {"domain", 'd', "d_x,d_y,d_z", 0, "gives the dimension of the domain (integers: d_x,d_y,d_z)"},
    {"patch", 10002, "o_x,o_y,o_z,l_x,l_y,l_z,lvl", 0, "indicate a patch of origin (floats: o_x,o_y,o_z), of length (floats: l_x,l_y,l_z) at level (integer: lvl)"},
    {0}};

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

static error_t parse_opt(int key, char* arg, struct argp_state* state) {
    parse_arg_t* arguments = reinterpret_cast<parse_arg_t*>(state->input);

    switch (key) {
        case 'd': {
            int*    length = arguments->length_;
            error_t err    = atoi_list(3, arg, length);
            m_log("domain length: %d %d %d", length[0], length[1], length[2]);
            return err;
        }
        case 10001: {
            bool*   period = arguments->period_;
            error_t err    = atob_list(3, arg, period);
            m_log("periodicity: %d %d %d", period[0], period[1], period[2]);
            return err;
        }
        case 10002: {
            real_t  data[7];
            error_t err   = atof_list(7, arg, data);
            lid_t   level = (lid_t)(data[6]);
            arguments->patch_.push_back(Patch(data, data + 3, level));
            m_log("patch: level %d starting (%f %f %f) of length (%f %f %f)", level, data[0], data[1], data[2], data[3], data[4], data[5]);
            return err;
        }
        default:
            return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

static struct argp argp = {options, parse_opt, 0, doc};

void ParseArgument(int argc, char** argv, parse_arg_t* arguments) {
    argp_parse(&argp, argc, argv, 0, 0, arguments);

    // string fname = "arguments.args";
    // int    rank;
    // MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // // rank 0 creates the folder if it doesn't exist
    // if (rank == 0) {
    //     struct stat st = {0};
    //     if (stat(fname.c_str(), &st) == -1) {
    //         mkdir(fname.c_str(), 0770);
    //     }
    // }
    // FILE* file =fopen(fname.c_str(),"w+");
    // if(file!= nullptr){
    //     fprintf()
    // }
}