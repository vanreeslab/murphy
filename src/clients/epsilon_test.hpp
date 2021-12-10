#ifndef SRC_CLIENTS_EPSILON_TEST_HPP
#define SRC_CLIENTS_EPSILON_TEST_HPP

#include "clients/testcase.hpp"
#include "core/pointers.hpp"
#include "core/types.hpp"

class EpsilonTest : public TestCase {

    bool break_me_ = false;

    level_t level_start_   = 4;
    level_t level_min_ = 0;
    level_t level_max_ = 10;

    real_t eps_start_ = 1.0;
    real_t delta_eps_ = 0.5;

   public:
    void InitParam(ParserArguments* param) override;
    void Run() override;
};

#endif