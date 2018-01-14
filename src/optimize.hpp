#pragma once
#include <functional>
#include "type.hpp"
#include "vector.hpp"

namespace scigg
{
// y, grad = f(x), f'(x)
typedef std::function<void(real_t& y, vector_mutable_view grad,
                           vector_const_view x)>
    diff1_func;

real_t quasi_newton_bfgs(const diff1_func& obj, vector_mutable_view x,
                         real_t& y, real_t first_step_len, real_t c1, real_t c2,
                         real_t tol, uint_t max_iter);
}
