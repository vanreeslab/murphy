#ifndef SRC_CLIENTS_EPSILON_TEST_HPP
#define SRC_CLIENTS_EPSILON_TEST_HPP

#include "clients/testcase.hpp"
#include "core/pointers.hpp"
#include "core/types.hpp"

class EpsilonTest : public TestCase {
    real_t  epsilon_start_ = 1.0;
    level_t level_start_   = 4;

   public:
    void InitParam(ParserArguments* param) override;
    void Run() override;
};

#endif