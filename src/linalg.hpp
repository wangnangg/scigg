#pragma once
#include "matrix.hpp"
namespace markovgg
{
// QR decomposition using Modified Gram-Schmidt
// A will be replaced by Q
void qr_decomp_mgs(matrix& AQ, matrix& R);

// QR decomposition using Householder Reflection
void qr_decomp_hr();

// QR decomposition using Givens Rotation
void qr_decomp_gr();
}
