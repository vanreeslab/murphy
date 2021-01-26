#ifndef SRC_CLIENTS_SIMPLEADVECTION_HPP
#define SRC_CLIENTS_SIMPLEADVECTION_HPP

#include "clients/testcase.hpp"
#include "core/pointers.hpp"
#include "grid/field.hpp"
#include "grid/grid.hpp"
#include "tools/parser.hpp"
#include "tools/prof.hpp"

class SimpleAdvection : public TestCase {
    m_ptr<Field> vel_;
    m_ptr<Field> scal_;
    m_ptr<Field> sol_;
    m_ptr<Grid>  grid_;
    m_ptr<Prof>  prof_;

    std::string folder_diag_ = "data";

   public:
    ~SimpleAdvection();
    void InitParam(ParserArguments* param) override;
    void Run() override;

   protected:
    void Diagnostics(const real_t time, const real_t dt, const lid_t iter);
};

#endif  // SRC_CLIENTS_SIMPLEADVECTION_HPP