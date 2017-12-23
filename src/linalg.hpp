#pragma once
#include "matrix.hpp"
namespace markovgg
{
// QR decomposition using Modified Gram-Schmidt
// A will be replaced by Q
void qr_decomp_mgs(matrix_mutable_view A, matrix_mutable_view R);

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
