#pragma once
#include <cassert>
#include <type.hpp>
#include <vector>
namespace markovgg
{
class vector_const_view
{
protected:
    const real_t *_const_data;
    size_t _inc;
    size_t _dim;

public:
    vector_const_view(const real_t *const_data, size_t dim, size_t inc)
        : _const_data(const_data), _inc(inc), _dim(dim)
    {
    }
    const real_t &operator[](size_t i) const
    {
        assert(i < _dim);
        return _const_data[i * _inc];
    }
    size_t dim() const { return _dim; }
    size_t inc() const { return _inc; }
};
class vector_mutable_view : public vector_const_view
{
protected:
    real_t *_mutable_data;

public:
    vector_mutable_view(real_t *mutable_data, size_t dim, size_t inc)
        : vector_const_view(mutable_data, dim, inc), _mutable_data(mutable_data)
    {
    }
    const real_t &operator[](size_t i) const
    {
        assert(i < dim());
        return _const_data[i * inc()];
    }
    real_t &operator[](size_t i)
    {
        assert(i < dim());
        return _mutable_data[i * inc()];
    }
};
class vector : public vector_mutable_view
{
public:
    vector(vector &&) = default;
    vector &operator=(vector &&) = default;
    vector(const vector &v) : vector_mutable_view(nullptr, v.dim(), v.inc())
    {
        _data = v._data;
        _const_data = &_data[0];
        _mutable_data = &_data[0];
    }
    vector &operator=(const vector &v)
    {
        _data = v._data;
        _const_data = &_data[0];
        _mutable_data = &_data[0];
        _inc = v._inc;
        _dim = v._dim;
    }

    vector(size_t dim, real_t val = 0.0)
        : vector_mutable_view(nullptr, dim, 1), _data(dim, val)
    {
        _const_data = &_data[0];
        _mutable_data = &_data[0];
    }

private:
    std::vector<real_t> _data;
};

bool operator==(const vector_const_view &x, const vector_const_view &y);
}
