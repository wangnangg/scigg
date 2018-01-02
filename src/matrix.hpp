#pragma once
#include <cassert>
#include <vector>
#include "type.hpp"
#include "vector.hpp"

namespace markovgg
{
class matrix_base
{
protected:
    size_t _m;
    size_t _n;
    size_t _ldim;
    size_t _row_step;
    size_t _col_step;
    bool _col_major;

    matrix_base(size_t m, size_t n, size_t ldim, bool col_major)
        : _m(m), _n(n), _ldim(ldim), _col_major(col_major)
    {
        if (col_major)
        {
            _row_step = 1;
            _col_step = ldim;
        }
        else
        {
            _row_step = ldim;
            _col_step = 1;
        }
    }

public:
    size_t ldim() const { return _ldim; }
    size_t m() const { return _m; }
    size_t n() const { return _n; }
    bool is_col_major() const { return _col_major; }
};

class matrix;
class matrix_mutable_view;
class matrix_const_view;

class matrix_const_view : public matrix_base
{
protected:
    const real_t* _const_data;

public:
    matrix_const_view(const real_t* data, size_t m, size_t n, size_t ldim,
                      bool col_major)
        : matrix_base(m, n, ldim, col_major), _const_data(data)
    {
    }
    matrix_const_view(const matrix_mutable_view& M);
    matrix_const_view(const matrix& M);
    const real_t& operator()(size_t i, size_t j) const
    {
        assert(i < m() && j < n());
        return _const_data[i * _row_step + j * _col_step];
    }
    vector_const_view row(size_t row_idx) const;
    vector_const_view col(size_t col_idx) const;
    matrix_const_view sub(size_t start_row, size_t start_col,
                          size_t end_row = 0, size_t end_col = 0) const;
    matrix_const_view transpose() const
    {
        return matrix_const_view(_const_data, n(), m(), ldim(), !_col_major);
    }
};

class matrix_mutable_view : public matrix_base
{
    real_t* _mutable_data;

public:
    matrix_mutable_view(real_t* mutable_data, size_t m, size_t n, size_t ldim,
                        bool col_major)
        : matrix_base(m, n, ldim, col_major), _mutable_data(mutable_data)
    {
    }
    matrix_mutable_view(matrix& M);
    real_t& operator()(size_t i, size_t j) const
    {
        assert(i < m() && j < n());
        return _mutable_data[i * _row_step + j * _col_step];
    }
    vector_mutable_view row(size_t row_idx) const;
    vector_mutable_view col(size_t col_idx) const;
    matrix_mutable_view sub(size_t start_row, size_t start_col,
                            size_t end_row = 0, size_t end_col = 0) const;
    matrix_mutable_view transpose() const
    {
        return matrix_mutable_view(_mutable_data, n(), m(), ldim(),
                                   !_col_major);
    }
};

class matrix : public matrix_base
{
    std::vector<real_t> _data;  // row major
public:
    matrix(size_t m, size_t n, real_t val = 0.0, bool col_major = false)
        : matrix_base(m, n, col_major ? m : n, col_major), _data(m * n, val)
    {
    }
    explicit matrix(matrix_const_view mat);
    const real_t& operator()(size_t i, size_t j) const
    {
        assert(i < m() && j < n());
        return _data[i * _row_step + j * _col_step];
    }
    real_t& operator()(size_t i, size_t j)
    {
        assert(i < m() && j < n());
        return _data[i * _row_step + j * _col_step];
    }
    vector_const_view row(size_t row_idx) const;
    vector_mutable_view row(size_t row_idx);
    vector_const_view col(size_t col_idx) const;
    vector_mutable_view col(size_t col_idx);
    matrix_const_view sub(size_t start_row, size_t start_col,
                          size_t end_row = 0, size_t end_col = 0) const;
    matrix_mutable_view sub(size_t start_row, size_t start_col,
                            size_t end_row = 0, size_t end_col = 0);
    matrix_const_view transpose() const
    {
        return matrix_const_view(&operator()(0, 0), n(), m(), ldim(),
                                 !_col_major);
    }
    matrix_mutable_view transpose()
    {
        return matrix_mutable_view(&operator()(0, 0), n(), m(), ldim(),
                                   !_col_major);
    }
};

bool operator==(matrix_const_view m1, matrix_const_view m2);

matrix create_matrix(size_t m, size_t n, const std::vector<double>& v);

matrix identity_matrix(size_t m);
}
