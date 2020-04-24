#ifndef SRC_STENCIL_HPP_
#define SRC_STENCIL_HPP_

#include "grid.hpp"
#include "murphy.hpp"
#include "operator.hpp"

class Stencil : public OperatorF2F {
    sid_t ida_   = 0;     //!< current source dimension
    bool  inner_ = true;  //!< application mode: true = inner, false = outer

    Grid* grid_;

   public:
    Stencil(Grid* grid);

    void ApplyOpF2F(const qid_t* qid, GridBlock* block, const Field* fid_src, Field* fid_trg) override;
    void operator()(Field* field_src, Field* field_trg);

   protected:
    virtual void ApplyOpDerivInner(const qid_t* qid, GridBlock* block, const Field* fid_src, Field* fid_trg) = 0;
    virtual void ApplyOpDerivOuter(const qid_t* qid, GridBlock* block, const Field* fid_src, Field* fid_trg) = 0;
};

#endif  // SRC_STENCIL_HPP_