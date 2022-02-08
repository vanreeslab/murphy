#ifndef SRC_CLIENTS_CONVERGENCE_WENO_HPP_
#define SRC_CLIENTS_CONVERGENCE_WENO_HPP_

#include "clients/testcase.hpp"

class ConvergenceWeno : public TestCase {
    bool    adapt_;
    bool    fix_weno_;
    level_t ilevel_;

    level_t level_min_;
    level_t level_max_;

    real_t eps_start_;
    real_t delta_eps_;

   public:
    ~ConvergenceWeno();
    void InitParam(ParserArguments* param) override;
    void Run() override;
};

#endif  // SRC_CLIENTS_SIMPLEADVECTION_HPP_