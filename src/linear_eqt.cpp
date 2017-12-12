#include "linear_eqt.hpp"
#include "debug_utils.hpp"
#include "matrix.hpp"
#include "matrix_oper.hpp"
#include <cassert>
namespace markovgg
{
real_t max_diff(const vec& v1, const vec& v2)
{
    assert(v1.dim() == v2.dim());
    real_t max = 0;
    for (size_t i = 0; i < v1.dim(); i++)
    {
        real_t diff = absolute(v1[i] - v2[i]);
        if (diff > max)
        {
            max = diff;
        }
    }
    return max;
}

real_t prec_check(const vec& x1, const vec& x2, const vec& b1, const vec& b2)
{
    real_t diff1 = max_diff(x1, x2);
    real_t diff2 = max_diff(b1, b2);
    return diff1 > diff2 ? diff1 : diff2;
}

real_t power_method(vec& x, const sqr_mat& A, iter_ctrl ctrl)
{
    set_vec_sum(x, 1.0);
    vec x_next(x.dim(), 0.0);
    vec b(x.dim(), 0.0);
    vec b_curr(x.dim(), 0.0);
    for (int_t ii = 0; ii < ctrl.max_iter / ctrl.check_interval; ii++)
    {
        for (int_t jj = 0; jj < ctrl.check_interval; jj++)
        {
            mul(x_next, x, A);
            set_vec_sum(x_next, 1.0);
            std::swap(x_next, x);
        }
        mul(b_curr, x, A);
        real_t prec = prec_check(x_next, x, x, b_curr);
        if (prec < ctrl.target_prec)
        {
            return prec;
        }
    }
    mul(b_curr, x, A);
    return prec_check(x_next, x, x, b_curr);
}

real_t sor_method(vec& x, const sqr_mat& A, const vec& b, real_t w,
                  iter_ctrl ctrl)
{
    assert(x.dim() == A.dim());
    assert(x.dim() == b.dim());
    if (near_zero(b, ctrl.target_prec))
    {
        fill(x, 0.0);
        return true;
    }
    vec x_next(x.dim(), 0.0);
    vec b_curr(x.dim(), 0.0);
    for (int_t ii = 0; ii < ctrl.max_iter / ctrl.check_interval; ii++)
    {
        for (int_t jj = 0; jj < ctrl.check_interval; jj++)
        {
            for (size_t i = 0; i < A.dim(); i++)
            {
                real_t diag = 0;
                real_t remain = b[i];
                for (const auto& e : A.col(i))
                {
                    if (e.idx < i)
                    {
                        remain -= e.val * x_next[e.idx];
                    }
                    else if (e.idx > i)
                    {
                        remain -= e.val * x[e.idx];
                    }
                    else
                    {
                        diag = e.val;
                    }
                }
                x_next[i] = (1 - w) * x[i] + w / diag * remain;
            }
            std::swap(x_next, x);
        }
        mul(b_curr, x, A);
        real_t prec = prec_check(x_next, x, b, b_curr);
        if (prec < ctrl.target_prec)
        {
            return prec;
        }
    }
    mul(b_curr, x, A);
    return prec_check(x_next, x, b, b_curr);
}

real_t sor_method(vec& x, const sqr_mat& A, const vec& b, real_t sum, real_t w,
                  iter_ctrl ctrl)
{
    assert(x.dim() == A.dim());
    assert(x.dim() == b.dim());
    vec x_next(x.dim(), 0.0);
    vec b_curr(x.dim(), 0.0);
    size_t last_col = A.dim() - 1;
    for (int_t ii = 0; ii < ctrl.max_iter / ctrl.check_interval; ii++)
    {
        for (int_t jj = 0; jj < ctrl.check_interval; jj++)
        {
            real_t last_remain = sum;
            for (size_t i = 0; i < last_col; i++)
            {
                real_t diag = 0;
                real_t remain = b[i];
                for (const auto& e : A.col(i))
                {
                    if (e.idx < i)
                    {
                        remain -= e.val * x_next[e.idx];
                    }
                    else if (e.idx > i)
                    {
                        remain -= e.val * x[e.idx];
                    }
                    else
                    {
                        diag = e.val;
                    }
                }
                x_next[i] = (1 - w) * x[i] + w / diag * remain;
                last_remain -= x_next[i];
            }
            x_next[last_col] = (1 - w) * x[last_col] + w * last_remain;
            std::swap(x_next, x);
        }
        mul(b_curr, x, A);
        real_t prec = prec_check(x_next, x, b, b_curr);
        if (prec < ctrl.target_prec)
        {
            return prec;
        }
    }
    mul(b_curr, x, A);
    return prec_check(x_next, x, b, b_curr);
}

}
