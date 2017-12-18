#pragma once
#include <cassert>
#include <vector>
#include "type.hpp"

namespace markovgg
{
class matrix
{
    std::vector<real_t> _data;  // row major
    size_t _m;                  // m by n
    size_t _n;

public:
    matrix(size_t m, size_t n, real_t val = 0.0)
        : _data(m * n, val), _m(m), _n(n)
    {
    }

    real_t &operator()(size_t i, size_t j)
    {
        assert(i < _m && j < _n);
        return _data[i * _n + j];
    }
    real_t operator()(size_t i, size_t j) const
    {
        assert(i < _m && j < _n);
        return _data[i * _n + j];
    }

    size_t m() const { return _m; }
    size_t n() const { return _n; }
    real_t *begin() { return &*_data.begin(); }
    real_t *end() { return &*_data.end(); }
    const real_t *begin() const { return &*_data.begin(); }
    const real_t *end() const { return &*_data.end(); }
};

bool operator==(const matrix &m1, const matrix &m2);
}
