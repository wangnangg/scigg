#include <cassert>
#include <sstream>
#include "blas.hpp"
#include "linalg.hpp"
#include "matvec_oper.hpp"
#include "spblas.hpp"
#include "splinalg.hpp"
#include "spmatrix_oper.hpp"
#include "vector.hpp"

namespace markovgg
{
void spdecomp_ilu(spmatrix_mutable_view A)
{
    assert(A.m() == A.n());
    assert(A.is_compressed_row());
}
}
