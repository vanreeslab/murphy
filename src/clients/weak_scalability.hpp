#ifndef SRC_CLIENTS_WEAKSCALABILITY_HPP_
#define SRC_CLIENTS_WEAKSCALABILITY_HPP_

#include "clients/testcase.hpp"
#include "grid/field.hpp"
#include "grid/grid.hpp"
#include "tools/parser.hpp"
#include "tools/prof.hpp"

class WeakScalability : public TestCase {
    Prof*       prof_;
    std::string folder_diag_ = "data";

   public:
    ~WeakScalability();
    void InitParam(ParserArguments* param) override;
    void Run() override;
};

#endif  // SRC_CLIENTS_WEAKSCALABILITY_HPP_