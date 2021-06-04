#ifndef SRC_TESTCASE_HPP_
#define SRC_TESTCASE_HPP_

#include "core/types.hpp"
#include "tools/parser.hpp"

class TestCase {
    // general testcase parameters
    bool dump_detail_ = false;

    iter_t iter_max_   = 0;
    iter_t iter_adapt_ = 0;
    iter_t iter_diag_  = 0;
    iter_t iter_dump_  = 0;

   public:
    bool dump_detail() const { return dump_detail_; }

    iter_t iter_max() const { return iter_max_; }
    iter_t iter_adapt() const { return iter_adapt_; }
    iter_t iter_diag() const { return iter_diag_; }
    iter_t iter_dump() const { return iter_dump_; }

    virtual ~TestCase(){};  // to force the call of the children class destructor

    virtual void InitParam(ParserArguments* param);
    virtual void Run() = 0;
};

#endif  // SRC_TESTCASE_HPP_