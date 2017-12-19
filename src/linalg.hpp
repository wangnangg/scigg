#pragma once
#include "matrix.hpp"
namespace markovgg
{
// QR decomposition using Modified Gram-Schmidt
// A will be replaced by Q
void QR_decomp_MGS(matrix& AQ, matrix& R);

// QR decomposition using Householder Reflection
void QR_decomp_HR();

// QR decomposition using Givens Rotation
void QR_decomp_GR();
}
