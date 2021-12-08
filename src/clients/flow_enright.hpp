#ifndef SRC_CLIENTS_FLOWENRIGHT_HPP_
#define SRC_CLIENTS_FLOWENRIGHT_HPP_

#include "clients/time_dependent.hpp"
#include "grid/field.hpp"

class FlowEnright : public TimeDependent {
    int  weno_     = 3;
    bool fix_weno_ = false;

    real_t cfl_ = 0.0;

    Field* scal_;
    Field* vel_;

    RKFunctor* advection_ = nullptr;
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