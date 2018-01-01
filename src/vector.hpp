#pragma once
#include <cassert>
#include <type.hpp>
#include <vector>
namespace markovgg
{
class vector_base
{
protected:
    size_t _inc;
    size_t _dim;
    vector_base(size_t inc, size_t dim) : _inc(inc), _dim(dim) {}

public:
    size_t dim() const { return _dim; }
    size_t inc() const { return _inc; }
};

class vector;
class vector_mutable_view;
class vector_const_view;

class vector_const_view : public vector_base
{
protected:
    const real_t *_const_data;

public:
    vector_const_view(const real_t *const_data, size_t dim, size_t inc)
        : vector_base(inc, dim), _const_data(const_data)
    {
    }
    vector_const_view(const vector &v);
    vector_const_view(const vector_mutable_view &v);
    const real_t &operator[](size_t i) const
    {
        assert(i < dim());
        return _const_data[i * inc()];
    }
    vector_const_view sub(size_t start, size_t end = 0) const
    {
        assert(end <= dim());
        if (end == 0)
        {
            end = dim();
        }
        return vector_const_view(&operator[](start), end - start, inc());
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
    vector_mutable_view(vector &v);
    real_t &operator[](size_t i) const
    {
        assert(i < dim());
        return _mutable_data[i * inc()];
    }
    vector_mutable_view sub(size_t start, size_t end = 0) const
    {
        assert(end <= dim());
        if (end == 0)
        {
            end = dim();
        }
        return vector_mutable_view(&operator[](start), end - start, inc());
    }
};

class vector : public vector_base
{
    std::vector<real_t> _data;

public:
    vector(size_t dim, real_t val = 0.0) : vector_base(1, dim), _data(dim, val)
    {
    }
    explicit vector(vector_const_view x);
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
    vector_mutable_view sub(size_t start, size_t end)
    {
        assert(end <= dim());
        if (end == 0)
        {
            end = dim();
        }
        return vector_mutable_view(&_data[0], end - start, inc());
    }
    vector_const_view sub(size_t start, size_t end) const
    {
        assert(end <= dim());
        if (end == 0)
        {
            end = dim();
        }
        return vector_const_view(&operator[](start), end - start, inc());
    }
};

bool operator==(vector_const_view x, vector_const_view y);

vector create_vector(size_t n, const std::vector<double> &v);
}
