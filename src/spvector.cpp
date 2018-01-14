#include "spvector.hpp"

namespace scigg
{
spvec_const_view spvec_const_view::masked(size_t start)
{
    size_t counter = nnz;
    for (size_t i = 0; i < nnz; i++)
    {
        if (idx[i] >= start)
        {
            break;
        }
        else
        {
            counter -= 1;
        }
    }
    return spvec_const_view(counter, idx + nnz - counter, val + nnz - counter);
}

spvec_mutable_view spvec_mutable_view::masked(size_t start)
{
    size_t counter = nnz;
    for (size_t i = 0; i < nnz; i++)
    {
        if (idx[i] >= start)
        {
            break;
        }
        else
        {
            counter -= 1;
        }
    }
    return spvec_mutable_view(counter, idx + nnz - counter,
                              val + nnz - counter);
}

const real_t* search_spvec_entry(spvec_const_view view, size_t idx)
{
    const size_t* begin = view.idx;
    const size_t* end = view.idx + view.nnz;
    const size_t* first = std::lower_bound(begin, end, idx);
    if (first != end && *first == idx)
    {
        return &view.val[first - begin];
    }
    else
    {
        return nullptr;
    }
}

real_t* search_spvec_entry(spvec_mutable_view view, size_t idx)
{
    const size_t* begin = view.idx;
    const size_t* end = view.idx + view.nnz;
    const size_t* first = std::lower_bound(begin, end, idx);
    if (first != end && *first == idx)
    {
        return &view.val[first - begin];
    }
    else
    {
        return nullptr;
    }
}
}
