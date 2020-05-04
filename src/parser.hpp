#ifndef SRC_PARSER_HPP_
#define SRC_PARSER_HPP_

#include <list>

#include "patch.hpp"

using std::list;

typedef struct parse_arg_t {
    bool        period_[3] = {false, false, false};  //!< the periodicity of the domain
    int         length_[3] = {1, 1, 1};              //!< the aspect ratio of the domain
    list<Patch> patch_;                              //!< list of imposed level regions for the initialization
} parse_arg_t;

void ParseArgument(int argc, char** argv, parse_arg_t* arguments);

#endif
