#ifndef SRC_CLIENTS_FLOWENRIGHT_HPP_
#define SRC_CLIENTS_FLOWENRIGHT_HPP_

#include "clients/time_dependent.hpp"
#include "grid/field.hpp"
#include "operator/advection.hpp"
#include "time/rkfunctor.hpp"

/**
 * @brief as the velocity is a function of time, we need to overwrite the value of the field before calling the rhs
 * 
 * To do so we simply reset the velocity before calling the RhsSet function
 * 
 */
class EnrightRhs: public RKFunctor{
    real_t time_final_;
    RKFunctor* advection_;
    Field* vel_;

    public:

    EnrightRhs() = delete;
    EnrightRhs(const real_t time, RKFunctor* advection,Field* vel);

    real_t cfl_rk3() const override {return advection_->cfl_rk3();}
    real_t rdiff_rk3() const override {return advection_->rdiff_rk3();}

    void RhsSet(Grid* grid, const real_t time, Field* field_u, Field* field_y) override; 
    void RhsAcc(Grid* grid, const real_t time, Field* field_u, Field* field_y) override;
};

class FlowEnright : public TimeDependent {
    int  weno_     = 3;
    bool fix_weno_ = false;

    real_t cfl_ = 0.0;

    Field* scal_;
    Field* vel_;

    RKFunctor* advection_ = nullptr;
    EnrightRhs* enright_rhs_ = nullptr;
    RK3_TVD*   rk3_       = nullptr;

    real_t time_accum_field_ = 0.0;
    real_t time_dump_field_  = 0.04;

    ~FlowEnright();

    std::string name() override { return "FlowEnright"; };



   protected:
    void Setup(ParserArguments* param) override;
    void DoTimeStep(real_t* time, real_t* dt) override;
    void Adapt(const real_t time, const real_t dt) override;
    void Diagnostics(const real_t time, const real_t dt, const lid_t iter, const real_t wtime) override;
};

#endif