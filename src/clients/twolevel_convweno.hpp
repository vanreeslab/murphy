#ifndef SRC_CLIENTS_TWOLEVEL_CONVWENO_HPP_
#define SRC_CLIENTS_TWOLEVEL_CONVWENO_HPP_

#include "clients/testcase.hpp"

class TwoLevelConvWeno : public TestCase {
    bool    adapt_;
    bool    fix_weno_;
    level_t ilevel_;

    level_t level_min_;
    level_t level_max_;

    real_t nu_;

   public:
    ~TwoLevelConvWeno();
    void InitParam(ParserArguments* param) override;
    void Run() override;
};

#endif  // SRC_CLIENTS_SIMPLEADVECTION_HPP_