#pragma once
#include "matrix.hpp"
namespace markovgg
{
// QR decomposition using Modified Gram-Schmidt
// A will be replaced by Q
void qr_decomp_mgs(matrix_mutable_view A, matrix_mutable_view R);

// find vector v, with v[0] = 1, so that (I - tau v v') w zeros w[1:n]. v[1:n]
// is stored in w[1:n] and w[0] = (Px)[0]
real_t find_householder_vector(vector_mutable_view w);

// w = w - \tau vv'w, v[0] is ignored and v[0] = 1 is assumed. v will only be
// changed temporarily.
void apply_householder_reflector(vector_mutable_view w, real_t tau,
                                 vector_mutable_view v);

// A = A - \tau vv'A, v[0] is ignored and v[0] = 1 is assumed. v will only be
// changed temporarily.
void apply_householder_reflector(matrix_mutable_view A, real_t tau,
                                 vector_mutable_view v);

// QR decomposition using householder reflection
// the upper triangular of A will be replaced by R
// the strict lower part will contain all the reflection vectors with v[0] = 1
// implicitly
void qr_decomp_hr(matrix_mutable_view A, vector_mutable_view tau_vec);

// unpack Q, R from QR and tau, where the strict lower part of V are reflection
// vectors from QR decomposition using householder reflection
void unpack_qr(matrix_mutable_view QR, vector_const_view tau_vec,
               matrix_mutable_view Q, matrix_mutable_view R);

// QR decomposition using Givens Rotation
void qr_decomp_gr();
}
