#pragma once
#include "matrix.hpp"
namespace markovgg
{
// QR decomposition using Modified Gram-Schmidt
// A will be replaced by Q
void decomp_qr_mgs(matrix_mutable_view A, matrix_mutable_view R);

// find vector v, with v[0] = 1, so that (I - tau v v') w zeros w[1:n]. v[1:n]
// is stored in w[1:n] and w[0] = (Px)[0]
real_t find_householder_vector(vector_mutable_view w);

// w = w - \tau vv'w, v[0] is ignored and v[0] = 1 is assumed. v will only be
// changed temporarily.
void apply_householder_reflector(vector_mutable_view w, real_t tau,
                                 vector_const_view v);

// A = A - \tau vv'A, v[0] is ignored and v[0] = 1 is assumed. v will only be
// changed temporarily.
void apply_householder_reflector(matrix_mutable_view A, real_t tau,
                                 vector_const_view v);

// QR decomposition using householder reflection
// the upper triangular of A will be replaced by R
// the strict lower part will contain all the reflection vectors with v[0] = 1
// implicitly
void decomp_qr_hr(matrix_mutable_view A, vector_mutable_view tau_vec);

// unpack Q, R from QR and tau, where the strict lower part of V are reflection
// vectors from QR decomposition using householder reflection
void unpack_qr(matrix_const_view QR, vector_const_view tau_vec,
               matrix_mutable_view Q, matrix_mutable_view R);

void qt_dot_vector(matrix_const_view QR, vector_const_view tau,
                   vector_mutable_view v);
void q_dot_vector(matrix_const_view QR, vector_const_view tau,
                  vector_mutable_view v);

// solve L x = b, b will be replaced with solution
void solve_lower_tri(matrix_const_view L, vector_mutable_view b);

// solve U x = b, b will be replaced with solution
void solve_upper_tri(matrix_const_view U, vector_mutable_view b);

// solve min||b - A x|| using QR decomposition, b will be replaced with solution
void least_square_qr(matrix_mutable_view A, vector_mutable_view b);

// PA = LU, A will be replaced by L and U, the rows of P will be permuted
void decomp_lu(matrix_mutable_view A, matrix_mutable_view P);

// A will be replaced by U
void unpack_lu(matrix_mutable_view A, matrix_mutable_view L);

// solve A * x = b using LU factorization. b will be replaced by x
void solve_lu(matrix_mutable_view A, vector_mutable_view b);

// solve A * x = b using QR factorization. b will be replaced by x
void solve_qr(matrix_mutable_view A, vector_mutable_view b);
}
