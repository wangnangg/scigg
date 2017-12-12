#pragma once

#include <exception>
#include <functional>
#include <utility>
#include "debug_utils.hpp"
#include "matrix.hpp"
#include "matrix_oper.hpp"
#include "type.hpp"

namespace markovgg
{
struct iter_ctrl
{
    real_t target_prec;
    int_t check_interval;
    int_t max_iter;
};

enum linear_solve_method
{
    SOR,
};

/** Solve x * A = x using power method. x is the initial value and cannot be 0.
 * The result will be normalized to 1.
 * reached precision is returned.
 */
real_t power_method(vec& x, const sqr_mat& A, iter_ctrl ctrl);

/** Solve x * A = b using sor method. A must be non-singular.
 * x is the initial guess and can be 0.
 * reached precision is returned.
 */
real_t sor_method(vec& x, const sqr_mat& A, const vec& b, real_t w,
                  iter_ctrl ctrl);

/** Solve x * A = b using sor method. A must be of rank n - 1.
 * x is the initial guess and can be 0.
 * sum is the sum of x.
 * reached precision is returned.
 */
real_t sor_method(vec& x, const sqr_mat& A, const vec& b, real_t sum, real_t w,
                  iter_ctrl ctrl);

real_t solve_linear_eq(vec& x, const sqr_mat& A, const vec& b, iter_ctrl ctrl,
                       linear_solve_method method);
real_t solve_linear_eq(vec& x, const sqr_mat& A, const vec& b, real_t sum,
                       iter_ctrl ctrl, linear_solve_method method);
}

