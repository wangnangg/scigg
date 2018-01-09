#pragma once

#include <exception>
#include <functional>
#include <utility>
#include "debug_utils.hpp"
#include "matrix.hpp"
#include "type.hpp"

namespace markovgg
{
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
                         real_t tol, uint_t check_interval);

real_t spsolve_restart_gmres_gms(spmatrix_const_view A, vector_mutable_view x,
                                 vector_const_view b,
                                 size_t kdim,  // dim of krylov space
                                 real_t tol, uint_t max_iter,
                                 uint_t check_interval);
}
