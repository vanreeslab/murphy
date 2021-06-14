#ifndef SRC_CLIENTS_CONVERGENCE_WENO_HPP_
#define SRC_CLIENTS_CONVERGENCE_WENO_HPP_

#include "clients/testcase.hpp"

class ConvergenceWeno : public TestCase {
    bool    adapt_  = false;
    level_t ilevel_ = 1;

   public:
    ~ConvergenceWeno();
    void InitParam(ParserArguments* param) override;
    void Run() override;
};

#endif  // SRC_CLIENTS_SIMPLEADVECTION_HPP_