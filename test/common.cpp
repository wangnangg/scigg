#include "common.hpp"
#include <iostream>
#include "matrix.hpp"

using namespace markovgg;
void print_sep() { std::cout << "-------------------" << std::endl; }

vec create_vec(const std::vector<real_t>& val)
{
    vec v(val.size(), 0.0);
    for (size_t i = 0; i < v.dim(); i++)
    {
        v[i] = val[i];
    }
    return v;
}

sqr_mat create_sqr_mat(size_t dim, const std::vector<real_t>& val)
{
    sqr_mat_crtor mat_crt(dim);
    assert(dim * dim == val.size());
    for (size_t i = 0; i < dim; i++)
    {
        for (size_t j = 0; j < dim; j++)
        {
            mat_crt.add_entry(i, j, val[i * dim + j]);
        }
    }
    return mat_crt.create();
}

matrix create_matrix(size_t row_dim, size_t col_dim,
                     const std::vector<real_t>& val)
{
    matrix_crtor mat_crt(row_dim, col_dim);
    assert(row_dim * col_dim == val.size());
    for (size_t i = 0; i < row_dim; i++)
    {
        for (size_t j = 0; j < col_dim; j++)
        {
            mat_crt.add_entry(i, j, val[i * col_dim + j]);
        }
    }
    return mat_crt.create();
}
