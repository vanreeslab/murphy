#include "parser.hpp"

#include <argp.h>
#include <stdlib.h>
#include <string.h>

#include "murphy.hpp"

static char doc[] = "MUPHY - a multiresolution multiphysics framework.";

static struct argp_option options[] = {
    {"periodic", 'p', "0/1,0/1,0/1", 0, "gives the periodicity in x,y,z (1 = periodic, 0 = not periodic) "},
    {"domain", 'd', "int,int,in", 0, "gives the dimension of the domain x,y,z (integers)"},
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

static error_t atof_list(const int length, char* arg, bool* list) {
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
        case 'p': {
            bool*   period = arguments->period;
            error_t err    = atob_list(3, arg, period);
            m_log("periodicity: %d %d %d\n", period[0], period[1], period[2]);
            return err;
        }
        case 'd': {
            int*    length = arguments->length;
            error_t err    = atoi_list(3, arg, length);
            m_log("domain length: %d %d %d\n", length[0], length[1], length[2]);
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
}