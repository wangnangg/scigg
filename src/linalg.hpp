#pragma once
#include "matrix.hpp"
namespace markovgg
{
// QR decomposition using Modified Gram-Schmidt
// A will be replaced by Q
void qr_decomp_mgs(matrix_mutable_view A, matrix_mutable_view R);

// QR decomposition using householder reflection
// A will be replaced by R
// V will contain all the reflection vectors (in lower triangular part)
void qr_decomp_hr(matrix_mutable_view A, matrix_mutable_view V);

// reconstruct Q from V, where V are reflection vectors from QR decomposition
// using householder reflection
void recover_q_from_v(matrix_const_view V, matrix_mutable_view Q);

// QR decomposition using Givens Rotation
void qr_decomp_gr();
}
