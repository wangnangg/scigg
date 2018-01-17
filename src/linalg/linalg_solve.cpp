#include <cassert>
#include "debug.hpp"
#include "linalg.hpp"
#include "matvec.hpp"

namespace scigg
{
// solve L x = b, b will be replaced with solution
void solve_lower_tri(matrix_const_view L, vector_mutable_view b)
{
    assert(L.m() == L.n());
    assert(L.m() == b.dim());
    b[0] /= L(0, 0);
    for (size_t i = 1; i < L.m(); i++)
    {
        b[i] -= dot(L.row(i).sub(0, i), b.sub(0, i));
        b[i] /= L(i, i);
    }
}

// solve L x = b, b will be replaced with solution. Diagonal is implicitly 1.
void solve_lower_tri_diag1(matrix_const_view L, vector_mutable_view b)
{
    assert(L.m() == L.n());
    assert(L.m() == b.dim());
    for (size_t i = 1; i < L.m(); i++)
    {
        b[i] -= dot(L.row(i).sub(0, i), b.sub(0, i));
    }
}

// solve U x = b, b will be replaced with solution
void solve_upper_tri(matrix_const_view U, vector_mutable_view b)
{
    assert(U.m() == U.n());
    assert(U.m() == b.dim());
    b[U.m() - 1] /= U(U.m() - 1, U.m() - 1);
    for (size_t j = 1; j < U.m(); j++)
    {
        size_t i = U.m() - j - 1;
        b[i] -= dot(U.row(i).sub(i + 1, U.m()), b.sub(i + 1, U.m()));
        b[i] /= U(i, i);
    }
}

void solve_lu(matrix_mutable_view A, vector_mutable_view b)
{
    assert(A.m() == A.n());
    size_t N = A.m();
    matrix_mutable_view b_(&b[0], N, 1, N, true);
    decomp_lup(A, b_);
    solve_lower_tri_diag1(A, b);
    solve_upper_tri(A, b);
}

void solve_qr(matrix_mutable_view A, vector_mutable_view b)
{
    assert(A.m() == A.n());
    size_t N = A.m();
    vector tau(N);
    decomp_qr_hr(A, tau);
    qt_dot_vector(A, tau, b);
    solve_upper_tri(A, b);
}
}
