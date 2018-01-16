#include "spmatvec_oper.hpp"
#include "debug_utils.hpp"
#include "spblas.hpp"
namespace scigg
{
spmatrix create_spmatrix(size_t m, size_t n, const std::vector<double>& v,
                         bool is_row_compressed)
{
    assert(m * n == v.size());
    spmatrix_creator M(m, n);
    for (size_t i = 0; i < m; i++)
    {
        for (size_t j = 0; j < n; j++)
        {
            auto val = v[i * n + j];
            if (val != 0)
            {
                M.add_entry(i, j, v[i * n + j]);
            }
        }
    }
    return M.create(is_row_compressed);
}

void dot(spmatrix_const_view A, vector_const_view x, vector_mutable_view y)
{
    // y = alpha * A * x + beta * y
    spblas_matrix_vector(1.0, A, x, 0.0, y);
}

matrix spmatrix2dense(spmatrix_const_view A)
{
    matrix M(A.m(), A.n(), 0.0);
    if (A.is_compressed_row())
    {
        for (size_t i = 0; i < A.ldim(); i++)
        {
            auto view = A[i];
            for (size_t j = 0; j < view.nnz; j++)
            {
                size_t idx = view.idx[j];
                real_t val = view.val[j];
                M(i, idx) = val;
            }
        }
    }
    else
    {
        for (size_t i = 0; i < A.ldim(); i++)
        {
            auto view = A[i];
            for (size_t j = 0; j < view.nnz; j++)
            {
                size_t idx = view.idx[j];
                real_t val = view.val[j];
                M(idx, i) = val;
            }
        }
    }
    return M;
}

real_t dot(spvec_const_view x, vector_const_view y)
{
    real_t prod = 0.0;
    for (size_t i = 0; i < x.nnz; i++)
    {
        size_t idx = x.idx[i];
        real_t val = x.val[i];
        prod += y[idx] * val;
    }
    return prod;
}

static size_t count_none_zero(spvec_const_view av, spvec_const_view bv)
{
    size_t ja = 0;
    size_t jb = 0;
    size_t j = 0;
    while (ja < av.nnz && jb < bv.nnz)
    {
        size_t ida = av.idx[ja];
        size_t idb = bv.idx[jb];
        if (ida < idb)
        {
            ja += 1;
        }
        else if (ida > idb)
        {
            jb += 1;
        }
        else
        {
            ja += 1;
            jb += 1;
        }
        j += 1;
    }
    j += av.nnz - ja;
    j += bv.nnz - jb;
    return j;
}

spmatrix add(spmatrix_const_view A, spmatrix_const_view B)
{
    assert(A.m() == B.m());
    assert(A.n() == B.n());
    assert(A.is_compressed_row() == B.is_compressed_row());
    size_t nnz = 0;
    std::vector<size_t> counter(A.ldim(), 0);
    for (size_t i = 0; i < A.ldim(); i++)
    {
        auto av = A[i];
        auto bv = B[i];
        size_t count = count_none_zero(av, bv);
        counter[i] = count;
        nnz += count;
    }
    std::vector<size_t> ptr(A.ldim() + 1, 0);
    std::vector<size_t> idx(nnz, 0);
    std::vector<real_t> val(nnz, 0);
    size_t pos = 0;
    ptr[0] = pos;
    for (size_t i = 0; i < A.ldim(); i++)
    {
        pos += counter[i];
        ptr[i + 1] = pos;
    }
    for (size_t i = 0; i < A.ldim(); i++)
    {
        auto av = A[i];
        auto bv = B[i];
        size_t ja = 0;
        size_t jb = 0;
        size_t j = 0;
        size_t* idxp = &idx[ptr[i]];
        real_t* valp = &val[ptr[i]];
        while (ja < av.nnz && jb < bv.nnz)
        {
            size_t ida = av.idx[ja];
            real_t vala = av.val[ja];
            size_t idb = bv.idx[jb];
            real_t valb = bv.val[jb];
            if (ida < idb)
            {
                idxp[j] = ida;
                valp[j] = vala;
                ja += 1;
            }
            else if (ida > idb)
            {
                idxp[j] = idb;
                valp[j] = valb;
                jb += 1;
            }
            else
            {
                idxp[j] = ida;
                valp[j] = vala + valb;
                ja += 1;
                jb += 1;
            }
            j += 1;
        }
        while (ja < av.nnz)
        {
            size_t idx = av.idx[ja];
            real_t val = av.val[ja];
            idxp[j] = idx;
            valp[j] = val;
            j += 1;
            ja += 1;
        }
        while (jb < bv.nnz)
        {
            size_t idx = bv.idx[jb];
            real_t val = bv.val[jb];
            idxp[j] = idx;
            valp[j] = val;
            j += 1;
            jb += 1;
        }
        assert(ptr[i] + j == ptr[i + 1]);
    }
    return spmatrix(A.m(), A.n(), A.is_compressed_row(), std::move(ptr),
                    std::move(idx), std::move(val));
}
}
