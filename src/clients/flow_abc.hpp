#ifndef SRC_CLIENTS_FLOW_ABC_HPP
#define SRC_CLIENTS_FLOW_ABC_HPP

#include "clients/testcase.hpp"
#include "core/pointers.hpp"
#include "grid/field.hpp"
#include "grid/grid.hpp"
#include "tools/parser.hpp"
#include "tools/prof.hpp"

class FlowABC : public TestCase {
    m_ptr<Field> vel_;
    m_ptr<Field> scal_;
    m_ptr<Grid>  grid_;
    m_ptr<Prof>  prof_;

    std::string folder_diag_ = "data";

   public:
    ~FlowABC();
    void InitParam(ParserArguments* param) override;
    void Run() override;

   protected:
    void Diagnostics(const real_t time, const real_t dt, const lid_t iter);
};

#endif  // SRC_CLIENTS_FLOW_ABC_HPP