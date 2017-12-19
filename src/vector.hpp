#pragma once
#include <cassert>
#include <type.hpp>
#include <vector>
namespace markovgg
{
class vector_base
{
    size_t _inc;
    size_t _dim;

protected:
    vector_base(size_t inc, size_t dim) : _inc(inc), _dim(dim) {}

public:
    size_t dim() const { return _dim; }
    size_t inc() const { return _inc; }
};

class vector : public vector_base
{
    std::vector<real_t> _data;

public:
    vector(size_t dim, real_t val = 0.0) : vector_base(1, dim), _data(dim, val)
    {
    }
    const real_t &operator[](size_t i) const
    {
        assert(i < dim());
        return _data[i * inc()];
    }
    real_t &operator[](size_t i)
    {
        assert(i < dim());
        return _data[i * inc()];
    }
};

class vector_mutable_view : public vector_base
{
protected:
    real_t *_mutable_data;

public:
    vector_mutable_view(real_t *mutable_data, size_t dim, size_t inc)
        : vector_base(inc, dim), _mutable_data(mutable_data)
    {
    }
    vector_mutable_view(vector &v)
        : vector_mutable_view(&v[0], v.dim(), v.inc())
    {
    }
    real_t &operator[](size_t i) const
    {
        assert(i < dim());
        return _mutable_data[i * inc()];
    }
};

class vector_const_view : public vector_base
{
protected:
    const real_t *_const_data;

public:
    vector_const_view(const real_t *const_data, size_t dim, size_t inc)
        : vector_base(inc, dim), _const_data(const_data)
    {
    }
    vector_const_view(const vector &v)
        : vector_const_view(&v[0], v.dim(), v.inc())
    {
    }
    vector_const_view(vector_mutable_view v)
        : vector_const_view(&v[0], v.dim(), v.inc())
    {
    }
    const real_t &operator[](size_t i) const
    {
        assert(i < dim());
        return _const_data[i * inc()];
    }
};

bool operator==(vector_const_view x, vector_const_view y);
}
