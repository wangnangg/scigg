#include "optimize.hpp"
#include <cassert>
#include "debug_utils.hpp"
#include "matrix.hpp"
#include "matvec_oper.hpp"
namespace markovgg
{
// on input:
//  x contains current solution
//  y = f(x)
//  grad = df(x)
// on output:
//  x is the next point
//  y = f(x)
//  grad = df(x)
void line_search_wolfe(const diff1_func& obj, real_t c1, real_t c2,
                       vector_const_view p, real_t& alpha,  //
                       vector_const_view x0, real_t y0,
                       vector_const_view grad0,  //
                       vector_mutable_view x1, real_t& y1,
                       vector_mutable_view grad1)
{
    real_t step_upper_bound = 0.0;
    do
    {
        // x1 = alpha * p + x0
        copy(x0, x1);
        blas_axpy(alpha, p, x1);
        obj(y1, grad1, x1);
        real_t deval = dot(grad0, p);
        if (y1 > y0 + c1 * alpha * deval)  // descent is not enough, shrink step
        {
            step_upper_bound = alpha;
            alpha /= 2.0;
        }
        else if (dot(grad1, p) <
                 c2 * deval)  // curvature not large enough, increase step
        {
            if (step_upper_bound == 0.0)  // haven't found the upper bound
            {
                alpha *= 2.0;
            }
            else
            {
                alpha = (alpha + step_upper_bound) / 2.0;
            }
        }
        else
        {  // wolfe conditions are satisfied
            break;
        }
    } while (true);
}

// Hk = Hk - p*s*(y^T*H) - p*(H*y)*s^T + p*p*s*(y^T*H*y)*s^T + p * s * s^T, s =
// x1 - x0, y = f'(x1) - f'(x0)
void update_hmatrix(matrix_mutable_view H, vector_const_view s,
                    vector_const_view y)
{
    vector v = dot(H, y);
    real_t p = 1.0 / dot(y, s);
    blas_rank1(-p, s, v, H);
    blas_rank1(-p, v, s, H);
    blas_rank1(p * p * dot(v, y) + p, s, s, H);
}

real_t quasi_newton_bfgs(const diff1_func& obj, vector_mutable_view x,
                         real_t& y, real_t first_step_len, real_t c1, real_t c2,
                         real_t tol, uint_t max_iter)
{
    assert(c2 > c1);
    assert(1.0 > c2);
    assert(c1 > 0.0);

    size_t N = x.dim();
    vector grad0(N);
    vector grad1(N);
    auto x0 = vector(x);
    real_t y0;
    real_t y1;
    vector p(N);
    vector dx(N);
    vector dg(N);

    // bootstrap
    obj(y0, grad0, x0);
    real_t grad_norm = norm2(grad0);
    vector x1 = x0;
    blas_axpy(-first_step_len / grad_norm, grad0, x1);
    obj(y1, grad1, x1);
    matrix H = identity_matrix(N);
    sub(x1, x0, dx);
    sub(grad1, grad0, dg);
    scale(dot(dg, dx) / dot(dg, dg), H);
    // start iteration
    for (size_t i = 0; i < max_iter; i++)
    {
        update_hmatrix(H, dx, dg);
        blas_matrix_vector(-1.0, H, grad1, 0.0, p);
        real_t alpha = 1.0;
        std::swap(x0, x1);
        std::swap(grad0, grad1);
        std::swap(y0, y1);
        line_search_wolfe(obj, c1, c2, p, alpha, x0, y0, grad0, x1, y1, grad1);
        grad_norm = norm2(grad1);
        if (grad_norm <= tol)
        {
            break;
        }
        sub(x1, x0, dx);
        sub(grad1, grad0, dg);
    }
    copy(x1, x);
    y = y1;
    return grad_norm;
}
}
