#ifndef SRC_CLIENTS_TIMEDEPENDENT_HPP_
#define SRC_CLIENTS_TIMEDEPENDENT_HPP_

#include "clients/testcase.hpp"
#include "grid/grid.hpp"
#include "tools/parser.hpp"
#include "tools/prof.hpp"
#include "time/rk3_tvd.hpp"

class TimeDependent : public TestCase {
    protected:
    // parameters
    bool no_adapt_ = false;
    bool optimal_tol_ = false;

    real_t tstart_ = 0.0;
    real_t tfinal_ = 1.0;

    level_t level_min_ = 0;
    level_t level_max_ = P8EST_QMAXLEVEL;

    std::string folder_diag_ = "data";

    // datastructures
    Grid* grid_ = nullptr;
    Prof* prof_ = nullptr;

   public:
    virtual ~TimeDependent();

    virtual std::string name() = 0;

    void InitParam(ParserArguments* param) override;
    void Run() override;

   protected:
    virtual void Setup(ParserArguments* param)                                                         = 0;
    virtual void DoTimeStep(real_t* time, real_t* dt)                                                  = 0;
    virtual void Adapt(const real_t time, const real_t dt)                                             = 0;
    virtual void Diagnostics(const real_t time, const real_t dt, const lid_t iter, const real_t wtime) = 0;
};

#endif  // SRC_CLIENTS_SIMPLEADVECTION_HPP_