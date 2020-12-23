#ifndef SRC_TESTCASE_HPP_
#define SRC_TESTCASE_HPP_

#include "parser.hpp"

class TestCase {
   public:
    virtual ~TestCase(){}; // to force the call of the children class destructor
    virtual void Run()                             = 0;
    virtual void InitParam(ParserArguments* param) = 0;
};

#endif  // SRC_TESTCASE_HPP_