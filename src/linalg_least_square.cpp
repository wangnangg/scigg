#include <cassert>
#include "blas.hpp"
#include "debug_utils.hpp"
#include "linalg.hpp"
#include "matvec_oper.hpp"

namespace markovgg
{
// A * x = b
void least_square_qr(matrix_mutable_view A, vector_mutable_view b)
{
    assert(A.m() >= A.n());
    vector tau(A.n());
    decomp_qr_hr(A, tau);
    qt_dot_vector(A, tau, b);
    solve_upper_tri(A.sub(0, 0, A.n(), A.n()), b.sub(0, A.n()));
}
}
