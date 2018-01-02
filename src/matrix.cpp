#include "matrix.hpp"

namespace markovgg
{
matrix_mutable_view::matrix_mutable_view(matrix& M)
    : matrix_mutable_view(&M(0, 0), M.m(), M.n(), M.ldim(), M.is_col_major())
{
}

matrix_const_view::matrix_const_view(const matrix& M)
    : matrix_const_view(&M(0, 0), M.m(), M.n(), M.ldim(), M.is_col_major())
{
}
matrix_const_view::matrix_const_view(const matrix_mutable_view& M)
    : matrix_const_view(&M(0, 0), M.m(), M.n(), M.ldim(), M.is_col_major())
{
}

matrix::matrix(matrix_const_view mat)
    : matrix_base(mat.m(), mat.n(), mat.is_col_major() ? mat.m() : mat.n(),
                  mat.is_col_major()),
      _data(mat.m() * mat.n())
{
    for (size_t i = 0; i < mat.m(); i++)
    {
        for (size_t j = 0; j < mat.n(); j++)
        {
            this->operator()(i, j) = mat(i, j);
        }
    }
}

vector_mutable_view matrix_mutable_view::row(size_t row_idx) const
{
    return vector_mutable_view(&operator()(row_idx, 0), n(), _col_step);
}

vector_mutable_view matrix_mutable_view::col(size_t col_idx) const
{
    return vector_mutable_view(&operator()(0, col_idx), m(), _row_step);
}

matrix_mutable_view matrix_mutable_view::sub(size_t start_row, size_t start_col,
                                             size_t end_row,
                                             size_t end_col) const
{
    assert(end_row <= m());
    assert(end_col <= n());
    if (end_row == 0)
    {
        end_row = m();
    }
    if (end_col == 0)
    {
        end_col = n();
    }
    size_t m = end_row - start_row;
    size_t n = end_col - start_col;
    return matrix_mutable_view(&operator()(start_row, start_col), m, n, ldim(),
                               is_col_major());
}

vector_const_view matrix_const_view::row(size_t row_idx) const
{
    return vector_const_view(&operator()(row_idx, 0), n(), _col_step);
}

vector_const_view matrix_const_view::col(size_t col_idx) const
{
    return vector_const_view(&operator()(0, col_idx), m(), _row_step);
}

matrix_const_view matrix_const_view::sub(size_t start_row, size_t start_col,
                                         size_t end_row, size_t end_col) const
{
    assert(end_row <= m());
    assert(end_col <= n());
    if (end_row == 0)
    {
        end_row = m();
    }
    if (end_col == 0)
    {
        end_col = n();
    }
    size_t m = end_row - start_row;
    size_t n = end_col - start_col;
    return matrix_const_view(&operator()(start_row, start_col), m, n, ldim(),
                             is_col_major());
}

vector_const_view matrix::row(size_t row_idx) const
{
    return vector_const_view(&operator()(row_idx, 0), n(), _col_step);
}

vector_const_view matrix::col(size_t col_idx) const
{
    return vector_const_view(&operator()(0, col_idx), m(), _row_step);
}

matrix_const_view matrix::sub(size_t start_row, size_t start_col,
                              size_t end_row, size_t end_col) const
{
    assert(end_row <= m());
    assert(end_col <= n());
    if (end_row == 0)
    {
        end_row = m();
    }
    if (end_col == 0)
    {
        end_col = n();
    }
    size_t m = end_row - start_row;
    size_t n = end_col - start_col;
    return matrix_const_view(&operator()(start_row, start_col), m, n, ldim(),
                             is_col_major());
}

vector_mutable_view matrix::row(size_t row_idx)
{
    return vector_mutable_view(&operator()(row_idx, 0), n(), _col_step);
}

vector_mutable_view matrix::col(size_t col_idx)
{
    return vector_mutable_view(&operator()(0, col_idx), m(), _row_step);
}

matrix_mutable_view matrix::sub(size_t start_row, size_t start_col,
                                size_t end_row, size_t end_col)
{
    assert(end_row <= m());
    assert(end_col <= n());
    if (end_row == 0)
    {
        end_row = m();
    }
    if (end_col == 0)
    {
        end_col = n();
    }
    size_t m = end_row - start_row;
    size_t n = end_col - start_col;
    return matrix_mutable_view(&operator()(start_row, start_col), m, n, ldim(),
                               is_col_major());
}

bool operator==(matrix_const_view m1, matrix_const_view m2)
{
    if (m1.m() != m2.m() || m1.n() != m2.n())
    {
        return false;
    }
    for (size_t i = 0; i < m1.m(); i++)
    {
        for (size_t j = 0; j < m1.n(); j++)
        {
            if (m1(i, j) != m2(i, j))
            {
                return false;
            }
        }
    }
    return true;
}
matrix create_matrix(size_t m, size_t n, const std::vector<double>& v)
{
    assert(m * n == v.size());
    matrix M(m, n);
    for (size_t i = 0; i < m; i++)
    {
        for (size_t j = 0; j < n; j++)
        {
            M(i, j) = v[i * n + j];
        }
    }
    return M;
}

matrix identity_matrix(size_t m)
{
    matrix M(m, m);
    for (size_t i = 0; i < m; i++)
    {
        M(i, i) = 1.0;
    }
    return M;
}
}
