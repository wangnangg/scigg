#pragma once
#include <functional>
#include "matvec.hpp"
#include "type.hpp"

namespace scigg
{
// y, grad = f(x), f'(x)
typedef std::function<real_t(vector_const_view x, vector_mutable_view grad)>
    diff1_obj_func;

real_t quasi_newton_bfgs(const diff1_obj_func& obj, vector_mutable_view x,
                         real_t& y, real_t first_step_len, real_t c1, real_t c2,
                         real_t tol, uint_t max_iter);
}
