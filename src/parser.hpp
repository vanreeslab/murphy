#ifndef SRC_PARSER_HPP_
#define SRC_PARSER_HPP_

#include <argp.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include <list>

#include "patch.hpp"

using std::list;

// typedef struct parse_arg_t {
//     int         n_repeat_  = 1;                      //!< the number of times we repeat the test
//     int         init_lvl_  = 1;                      //!< the level at which the grid is initialized
//     bool        period_[3] = {false, false, false};  //!< the periodicity of the domain
//     int         length_[3] = {1, 1, 1};              //!< the aspect ratio of the domain
//     list<Patch> patch_;                              //!< list of imposed level regions for the initialization
//     int         eta_smooting_ = 5;                   //!< the number of smooting steps in the MG
// } parse_arg_t;

// void ParseArgument(int argc, char** argv, parse_arg_t* arguments);

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

#endif
