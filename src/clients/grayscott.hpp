#ifndef SRC_CLIENTS_GRAYSCOTT_HPP
#define SRC_CLIENTS_GRAYSCOTT_HPP

#include "clients/testcase.hpp"
#include "core/pointers.hpp"
#include "grid/field.hpp"
#include "grid/grid.hpp"
#include "tools/parser.hpp"
#include "tools/prof.hpp"

class GrayScott : public TestCase {
    

    Field* lvlset_  = nullptr;
    Field* vel_ = nullptr;
    Field* u_ = nullptr;
    Field* v_ = nullptr;
    Grid*  grid_ = nullptr;
    Prof*  prof_ = nullptr;

    std::string folder_diag_ = "data";

   public:
    ~FlowABC();
    void InitParam(ParserArguments* param) override;
    void Run() override;

   protected:
    void Diagnostics(const real_t time, const real_t dt, const lid_t iter);
};

#endif  // SRC_CLIENTS_FLOW_ABC_HPP