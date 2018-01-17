#pragma once
#include <cassert>
#include <vector>
#include "type.hpp"

namespace scigg
{
struct spvec_mutable_view
{
    size_t nnz;
    size_t* idx;
    real_t* val;
    spvec_mutable_view(size_t nnz, size_t* idx, real_t* val)
        : nnz(nnz), idx(idx), val(val)
    {
    }
    spvec_mutable_view masked(size_t start);
};

struct spvec_const_view
{
    size_t nnz;
    const size_t* idx;
    const real_t* val;
    spvec_const_view(size_t nnz, const size_t* idx, const real_t* val)
        : nnz(nnz), idx(idx), val(val)
    {
    }
    spvec_const_view masked(size_t start);
    spvec_const_view(const spvec_mutable_view& view)
        : nnz(view.nnz), idx(view.idx), val(view.val)
    {
    }
};

const real_t* search_spvec_entry(spvec_const_view view, size_t idx);
real_t* search_spvec_entry(spvec_mutable_view view, size_t idx);
}
