#include "autodiff.hpp"

namespace scigg
{
real_t diff1_sum::eval(vector_const_view x, vector_mutable_view grad) const
{
    assert(x.dim() == grad.dim());
    vector child_grad(grad.dim());
    fill(0.0, grad);
    real_t val = 0.0;
    for (const auto& c : _children)
    {
        val += c->eval(x, child_grad);
        grad += child_grad;
    }
    return val;
}

real_t diff1_var::eval(vector_const_view x, vector_mutable_view grad) const
{
    assert(x.dim() == grad.dim());
    assert(_var_idx < x.dim());
    fill(0.0, grad);
    grad[_var_idx] = 1.0;
    return x[_var_idx];
}

real_t diff1_prod::eval(vector_const_view x, vector_mutable_view grad) const
{
    assert(x.dim() == grad.dim());
    real_t val = 1.0;
    matrix child_grad(grad.dim(), _children.size());
    vector child_val(_children.size());
    vector amp(child_val.dim());
    for (size_t i = 0; i < _children.size(); i++)
    {
        child_val[i] = _children[i]->eval(x, child_grad.col(i));
        val *= child_val[i];
    }
    if (val > 0)  // fast path
    {
        for (size_t i = 0; i < child_val.dim(); i++)
        {
            amp[i] = val / child_val[i];
        }
    }
    else
    {
        for (size_t i = 0; i < child_val.dim(); i++)
        {
            amp[i] = 1.0;
            for (size_t j = 0; j < child_val.dim(); j++)
            {
                if (i == j) continue;
                amp[i] *= child_val[j];
            }
        }
    }
    dot(child_grad, amp, grad);
    return val;
}

real_t diff1_con::eval(vector_const_view x, vector_mutable_view grad) const
{
    assert(x.dim() == grad.dim());
    fill(0.0, grad);
    return _val;
}
}
