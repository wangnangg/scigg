#pragma once
#include <cassert>
#include <vector>
#include "type.hpp"
#include "vector.hpp"

namespace markovgg
{
class matrix_const_view
{
protected:
    const real_t* _const_data;
    size_t _m;
    size_t _n;
    size_t _ldim;

public:
    matrix_const_view(const real_t* const_data, size_t m, size_t n, size_t ldim)
        : _m(m), _n(n), _ldim(ldim), _const_data(const_data)
    {
    }
    const real_t& operator()(size_t i, size_t j) const
    {
        assert(i < _m && j < _n);
        return _const_data[i * _ldim + j];
    }
    size_t ldim() const { return _ldim; }
    size_t m() const { return _m; }
    size_t n() const { return _n; }
};

class matrix_mutable_view : public matrix_const_view
{
protected:
    real_t* _mutable_data;

public:
    matrix_mutable_view(real_t* mutable_data, size_t m, size_t n, size_t ldim)
        : matrix_const_view(mutable_data, m, n, ldim),
          _mutable_data(mutable_data)
    {
    }
    const real_t& operator()(size_t i, size_t j) const
    {
        assert(i < m() && j < n());
        return _const_data[i * ldim() + j];
    }
    real_t& operator()(size_t i, size_t j)
    {
        assert(i < m() && j < n());
        return _mutable_data[i * ldim() + j];
    }
};

class matrix : public matrix_mutable_view
{
    std::vector<real_t> _data;  // row major
public:
    matrix(size_t m, size_t n, real_t val = 0.0)
        : matrix_mutable_view(nullptr, m, n, n), _data(m * n, val)
    {
        _const_data = &_data[0];
        _mutable_data = &_data[0];
    }
    matrix(matrix&&) = default;
    matrix& operator=(matrix&&) = default;
    matrix(const matrix& M)
        : matrix_mutable_view(nullptr, M.m(), M.n(), M.ldim()), _data(0)
    {
        _data = M._data;
        _const_data = &_data[0];
        _mutable_data = &_data[0];
    }
    matrix& operator=(const matrix& M)
    {
        _data = M._data;
        _const_data = &_data[0];
        _mutable_data = &_data[0];
        _m = M._m;
        _n = M._n;
        _ldim = M._ldim;
    }
};

bool operator==(const matrix_const_view& m1, const matrix_const_view& m2);
}
