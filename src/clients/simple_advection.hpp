#ifndef SRC_CLIENTS_SIMPLEADVECTION_HPP_
#define SRC_CLIENTS_SIMPLEADVECTION_HPP_

#include "clients/testcase.hpp"
#include "grid/field.hpp"
#include "grid/grid.hpp"
#include "tools/parser.hpp"
#include "tools/prof.hpp"

class SimpleAdvection : public TestCase {
    bool no_adapt_;
    bool grid_on_sol_;
    int  weno_;
    bool fix_weno_;
    bool refine_only_ = false;
    bool optimal_tol_ = false;

    real_t tstart_ = 0.0;
    real_t tfinal_ = 0.0;

    real_t t_deterr_       = 0.25;
    real_t t_deterr_accum_ = 0.0;

    real_t cfl_ = 0.0;
    real_t nu_  = 0.0;

    level_t level_min_ = 0;
    level_t level_max_ = P8EST_QMAXLEVEL;

    Field* vel_;
    Field* scal_;
    Grid*  grid_;
    Prof*  prof_;

    std::string folder_diag_ = "data";

   public:
    ~SimpleAdvection();
    void InitParam(ParserArguments* param) override;
    void Run() override;

   protected:
    void Diagnostics(const real_t time, const real_t dt, const lid_t iter,const real_t wtime);
    void GridDetErr(const real_t time, const real_t dt, const lid_t iter,const real_t wtime);
};

#endif  // SRC_CLIENTS_SIMPLEADVECTION_HPP_