#include "linalg.hpp"
#include <cassert>
#include "blas.hpp"
#include "matvec_oper.hpp"
namespace markovgg
{
void QR_decomp_MGS(matrix& AQ, matrix& R)
{
    assert(AQ.n() <= AQ.m());
    assert(R.n() == R.m());
    assert(AQ.n() == R.m());
    for (size_t i = 0; i < AQ.n(); i++)
    {
        auto qi = col_mutable_view(AQ, i);
        real_t rii = blas_norm2(qi);
        blas_scale(1.0 / rii, qi);
        R(i, i) = rii;
        for (size_t j = i + 1; j < AQ.n(); j++)
        {
            auto vj = col_mutable_view(AQ, j);
            real_t rij = blas_dot(qi, vj);
            // vj = vj - rjj * qi
            blas_axpy(-rij, qi, vj);
            R(i, j) = rij;
        }
    }
}
}
