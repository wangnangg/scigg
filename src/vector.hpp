#pragma once
#include <type.hpp>
#include <vector>

namespace markovgg
{
class vector
{
public:
    vector(vector &&) = default;
    vector &operator=(vector &&) = default;
    vector(const vector &) = default;
    vector &operator=(const vector &) = default;
    vector(size_t dim, real_t val = 0.0) : _data(dim, val) {}

    real_t operator[](size_t i) const { return _data[i]; }
    real_t &operator[](size_t i) { return _data[i]; }
    size_t dim() const { return _data.size(); }

    const real_t *begin() const { return &*_data.begin(); }
    const real_t *end() const { return &*_data.end(); }
    real_t *begin() { return &*_data.begin(); }
    real_t *end() { return &*_data.end(); }

private:
    std::vector<real_t> _data;
};

bool operator==(const vector &x, const vector &y);
}
