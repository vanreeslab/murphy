#include "testcase.hpp"

void TestCase::InitParam(ParserArguments* param) {
    //-------------------------------------------------------------------------
    dump_detail_ = param->dump_detail;
    iter_max_    = param->iter_max;
    iter_diag_   = param->iter_diag;
    iter_adapt_  = param->iter_adapt;
    iter_dump_   = param->iter_dump;
    //-------------------------------------------------------------------------
    m_log("testcase initialized with %d iter_max");
}