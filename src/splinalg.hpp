#pragma once

#include <exception>
#include <functional>
#include <utility>
#include "debug_utils.hpp"
#include "matrix.hpp"
#include "type.hpp"

namespace scigg
{
// implicit form of M^-1 solve for M x = b, x will be stored in b
typedef std::function<void(vector_mutable_view b)> pre_condition;

// incomplete LU factorization. L is stored in the strict lower part with
// diagonals implicitly equal 1.
void spdecomp_ilu(spmatrix_mutable_view A);

/** Solve A * x = x using power method. x is the initial value and cannot be 0.
 * The result will be normalized to 1.
 * reached precision is returned.
 */
real_t eigen_power_method(spmatrix_const_view A, vector_mutable_view x,
                          real_t tol, uint_t max_iter, uint_t check_interval);

/** Solve A * x = b using sor method. A must be non-singular.
 * x is the initial guess and can be 0.
 * reached precision is returned.
 */
real_t spsolve_sor_method(spmatrix_const_view A, vector_mutable_view x,
                          vector_const_view b, real_t w, real_t tol,
                          uint_t max_iter, uint_t check_interval);

real_t spsolve_gmres_gms(spmatrix_const_view A, vector_mutable_view x,
                         vector_const_view b,
                         size_t kdim,  // dim of krylov space
                         real_t tol, uint_t check_interval,
                         pre_condition Msolve = nullptr);

real_t spsolve_restart_gmres_gms(spmatrix_const_view A, vector_mutable_view x,
                                 vector_const_view b,
                                 size_t kdim,  // dim of krylov space
                                 real_t tol, uint_t max_iter,
                                 uint_t check_interval,
                                 pre_condition Msolve = nullptr);

// U * x = b; U is a upper triangular sparse matrix.
void spsolve_upper_tri(spmatrix_const_view U, vector_mutable_view b);
// L * x = b; L is a lower triangular sparse matrix.
void spsolve_lower_tri(spmatrix_const_view L, vector_mutable_view b);
void spsolve_lower_tri_diag1(spmatrix_const_view L, vector_mutable_view b);

// solve ilu: L * U * x = b
void spsolve_ilu(spmatrix_const_view iLU, vector_mutable_view b);
}
