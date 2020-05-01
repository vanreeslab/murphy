#ifndef SRC_PARSER_HPP_
#define SRC_PARSER_HPP_

typedef struct parse_arg_t {
    bool period[3] = {false,false,false}; //!< the periodicity of the domain
    int  length[3] = {1,1,1}; //!< the aspect ratio of the domain
}parse_arg_t;


void ParseArgument(int argc, char** argv,parse_arg_t* arguments) ;


#endif
