#ifndef SRC_CLIENTS_SIMPLEADVECTION_HPP_
#define SRC_CLIENTS_SIMPLEADVECTION_HPP_

#include "clients/testcase.hpp"
#include "core/pointers.hpp"
#include "grid/field.hpp"
#include "grid/grid.hpp"
#include "tools/parser.hpp"
#include "tools/prof.hpp"

class SimpleAdvection : public TestCase {
    bool no_adapt_;
    bool grid_on_sol_;
    int  weno_;

    real_t tstart_ = 0.0;
    real_t tfinal_ = 0.0;

    real_t cfl_ = 0.0;

    level_t level_min_ = 0;
    level_t level_max_ = P8EST_QMAXLEVEL;

    m_ptr<Field> vel_;
    m_ptr<Field> scal_;
    m_ptr<Field> sol_;
    m_ptr<Grid>  grid_;
    m_ptr<Prof>  prof_;

    m_ptr<SetScalarRing> ring_;
    m_ptr<SetPolynom>    vel_field_;

    std::string folder_diag_ = "data";

   public:
    ~SimpleAdvection();
    void InitParam(ParserArguments* param) override;
    void Run() override;

   protected:
    void Diagnostics(const real_t time, const real_t dt, const lid_t iter,const real_t wtime);
};

#endif  // SRC_CLIENTS_SIMPLEADVECTION_HPP_