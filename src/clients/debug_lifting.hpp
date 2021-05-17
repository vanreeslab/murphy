#ifndef SRC_CLIENTS_DEBUGLIFTING_HPP
#define SRC_CLIENTS_DEBUGLIFTING_HPP

#include "clients/testcase.hpp"
#include "core/pointers.hpp"
#include "core/types.hpp"

class DebugLifting : public TestCase {
   public:
    void InitParam(ParserArguments* param) override;
    void Run() override;
};

#endif