#ifndef SRC_PARSER_HPP_
#define SRC_PARSER_HPP_

#include <list>

#include "patch.hpp"

using std::list;

typedef struct parse_arg_t {
    int         n_repeat_  = 1; //!< the number of times we repeat the test
    int         init_lvl_  = 1; //!< the level at which the grid is initialized
    bool        period_[3] = {false, false, false};  //!< the periodicity of the domain
    int         length_[3] = {1, 1, 1};              //!< the aspect ratio of the domain
    list<Patch> patch_;                              //!< list of imposed level regions for the initialization
} parse_arg_t;

void ParseArgument(int argc, char** argv, parse_arg_t* arguments);

#endif