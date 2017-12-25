#pragma once
#include <cassert>
#include <vector>
#include "type.hpp"
#include "vector.hpp"

namespace markovgg
{
class matrix_base
{
    size_t _m;
    size_t _n;
    size_t _ldim;

protected:
    matrix_base(size_t m, size_t n, size_t ldim) : _m(m), _n(n), _ldim(ldim) {}

public:
    size_t ldim() const { return _ldim; }
    size_t m() const { return _m; }
    size_t n() const { return _n; }
};

class matrix;
class matrix_mutable_view;
class matrix_const_view;

class matrix_mutable_view : public matrix_base
{
protected:
    real_t* _mutable_data;

public:
    matrix_mutable_view(real_t* mutable_data, size_t m, size_t n, size_t ldim)
        : matrix_base(m, n, ldim), _mutable_data(mutable_data)
    {
    }
    matrix_mutable_view(matrix& M);
    real_t& operator()(size_t i, size_t j) const
    {
        assert(i < m() && j < n());
        return _mutable_data[i * ldim() + j];
    }
    vector_mutable_view row(size_t row_idx) const;
    vector_mutable_view col(size_t col_idx) const;
    matrix_mutable_view sub(size_t start_row, size_t start_col,
                            size_t end_row = 0, size_t end_col = 0) const;
};

class matrix_const_view : public matrix_base
{
protected:
    const real_t* _const_data;

public:
    matrix_const_view(const real_t* data, size_t m, size_t n, size_t ldim)
        : matrix_base(m, n, ldim), _const_data(data)
    {
    }
    matrix_const_view(matrix_mutable_view M);
    matrix_const_view(const matrix& M);
    const real_t& operator()(size_t i, size_t j) const
    {
        assert(i < m() && j < n());
        return _const_data[i * ldim() + j];
    }
    vector_const_view row(size_t row_idx) const;
    vector_const_view col(size_t col_idx) const;
    matrix_const_view sub(size_t start_row, size_t start_col,
                          size_t end_row = 0, size_t end_col = 0) const;
};

class matrix : public matrix_base
{
    std::vector<real_t> _data;  // row major
public:
    matrix(size_t m, size_t n, real_t val = 0.0)
        : matrix_base(m, n, n), _data(m * n, val)
    {
    }
    const real_t& operator()(size_t i, size_t j) const
    {
        assert(i < m() && j < n());
        return _data[i * ldim() + j];
    }
    real_t& operator()(size_t i, size_t j)
    {
        assert(i < m() && j < n());
        return _data[i * ldim() + j];
    }
    vector_const_view row(size_t row_idx) const;
    vector_mutable_view row(size_t row_idx);
    vector_const_view col(size_t col_idx) const;
    vector_mutable_view col(size_t col_idx);
    matrix_const_view sub(size_t start_row, size_t start_col,
                          size_t end_row = 0, size_t end_col = 0) const;
    matrix_mutable_view sub(size_t start_row, size_t start_col,
                            size_t end_row = 0, size_t end_col = 0);
};

bool operator==(matrix_const_view m1, matrix_const_view m2);

matrix create_matrix(size_t m, size_t n, const std::vector<double>& v);
}
