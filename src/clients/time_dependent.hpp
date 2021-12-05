#ifndef SRC_CLIENTS_TIMEDEPENDENT_HPP_
#define SRC_CLIENTS_TIMEDEPENDENT_HPP_

#include "clients/testcase.hpp"
#include "grid/field.hpp"
#include "grid/grid.hpp"
#include "tools/parser.hpp"
#include "tools/prof.hpp"

class TimeDependent : public TestCase {
    // parameters
    bool no_adapt_ = false;
    bool fix_weno_ = true;
    bool refine_only_ = false;
    bool optimal_tol_ = false;

    int  weno_ = 3;
    real_t tstart_ = 0.0;
    real_t tfinal_ = 0.0;

    level_t level_min_ = 0;
    level_t level_max_ = P8EST_QMAXLEVEL;

    // diagnostics time
    real_t t_deterr_       = 0.25;
    real_t t_deterr_accum_ = 0.0;

    std::string name_ = "time_dep_default";
    std::string folder_diag_ = "data";


    // datastructures
    Grid*  grid_;
    Prof*  prof_;

    

   public:
    ~TimeDependent();
    void InitParam(ParserArguments* param) override;
    void Run() override;

   protected:
    virtual void Setup(ParserArguments* param) = 0;
    virtual void InitCondition() = 0;
    void Diagnostics(const real_t time, const real_t dt, const lid_t iter,const real_t wtime);
};

#endif  // SRC_CLIENTS_SIMPLEADVECTION_HPP_