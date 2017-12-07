#pragma once
#include <cassert>
#include <string>
#include <vector>
#include "type.hpp"
namespace markovgg
{
class vec
{
public:
    vec(vec &&) = default;
    vec &operator=(vec &&) = default;
    vec(const vec &) = default;
    vec &operator=(const vec &) = default;
    vec(size_t dim, real_t val) : _val(dim, val) {}

    real_t operator[](size_t i) const { return _val[i]; }
    real_t &operator[](size_t i) { return _val[i]; }
    size_t dim() const { return _val.size(); }

    const real_t *begin() const { return &*_val.begin(); }
    const real_t *end() const { return &*_val.end(); }
    real_t *begin() { return &*_val.begin(); }
    real_t *end() { return &*_val.end(); }

private:
    std::vector<real_t> _val;
    friend bool operator==(const vec &x, const vec &y);
};

bool operator==(const vec &x, const vec &y);

struct sps_entry
{
    size_t idx;
    real_t val;
};

class sps_vec_view
{
public:
    sps_vec_view(const sps_entry *start, const sps_entry *end)
        : _start(start), _end(end)
    {
    }
    const sps_entry *begin() const { return _start; }
    const sps_entry *end() const { return _end; }
    size_t entry_count() const { return _end - _start; }

private:
    const sps_entry *_start;
    const sps_entry *_end;
};

class sps_vec_mut_view
{
public:
    sps_vec_mut_view(sps_entry *start, sps_entry *end)
        : _start(start), _end(end)
    {
    }
    sps_entry *begin() const { return _start; }
    sps_entry *end() const { return _end; }
    size_t entry_count() const { return _end - _start; }

private:
    sps_entry *_start;
    sps_entry *_end;
};

class matrix
{
protected:
    size_t _row_dim;
    std::vector<size_t> _idx_list;
    std::vector<sps_entry> _entry_list;

public:
    matrix(matrix &&) = default;
    matrix(const matrix &) = default;
    matrix &operator=(const matrix &) = default;
    matrix &operator=(matrix &&) = default;
    matrix(size_t row_dim, std::vector<size_t> idx_list,
           std::vector<sps_entry> entry_list)
        : _row_dim(row_dim),
          _idx_list(std::move(idx_list)),
          _entry_list(std::move(entry_list))
    {
    }

    size_t row_dim() const { return _row_dim; }
    size_t col_dim() const { return _idx_list.size() - 1; }
    size_t entry_count() const { return _entry_list.size(); }
    sps_vec_view col(size_t idx) const
    {
        size_t start_idx = _idx_list[idx];
        size_t size = _idx_list[idx + 1] - start_idx;
        const sps_entry *start_addr = &_entry_list[start_idx];
        const sps_entry *end_addr = start_addr + size;
        return sps_vec_view(start_addr, end_addr);
    }
    sps_vec_mut_view col(size_t idx)
    {
        size_t start_idx = _idx_list[idx];
        size_t size = _idx_list[idx + 1] - start_idx;
        sps_entry *start_addr = &_entry_list[start_idx];
        sps_entry *end_addr = start_addr + size;
        return sps_vec_mut_view(start_addr, end_addr);
    }
    virtual real_t operator()(size_t row_idx, size_t col_idx) const;
};

class sqr_mat : public matrix
{
public:
    sqr_mat(sqr_mat &&) = default;
    sqr_mat(const sqr_mat &) = default;
    sqr_mat &operator=(const sqr_mat &) = default;
    sqr_mat &operator=(sqr_mat &&) = default;
    size_t dim() const { return _row_dim; }
    size_t entry_count() const { return _entry_list.size(); }

    sqr_mat(size_t dim, std::vector<size_t> idx_list,
            std::vector<sps_entry> entry_list)
        : matrix(dim, std::move(idx_list), std::move(entry_list))
    {
    }

    sps_vec_view col(size_t idx) const
    {
        size_t start_idx = _idx_list[idx];
        size_t size = _idx_list[idx + 1] - start_idx;
        const sps_entry *start_addr = &_entry_list[start_idx];
        const sps_entry *end_addr = start_addr + size;
        return sps_vec_view(start_addr, end_addr);
    }
    sps_vec_mut_view col(size_t idx)
    {
        size_t start_idx = _idx_list[idx];
        size_t size = _idx_list[idx + 1] - start_idx;
        sps_entry *start_addr = &_entry_list[start_idx];
        sps_entry *end_addr = start_addr + size;
        return sps_vec_mut_view(start_addr, end_addr);
    }
    real_t diag(size_t idx) const
    {
        size_t start_idx = _idx_list[idx];
        return _entry_list[start_idx].val;
    }
    real_t &diag(size_t idx)
    {
        size_t start_idx = _idx_list[idx];
        return _entry_list[start_idx].val;
    }
    real_t operator()(size_t row_idx, size_t col_idx) const override;
};

struct mat_entry
{
    size_t row_idx;
    size_t col_idx;
    real_t val;
};

class matrix_crtor
{
protected:
    size_t _row_dim;
    size_t _col_dim;
    std::vector<mat_entry> _entry_list;

public:
    matrix_crtor(matrix_crtor &&) = default;
    matrix_crtor(size_t row_dim, size_t col_dim)
        : _row_dim(row_dim), _col_dim(col_dim)
    {
    }
    void add_entry(size_t row_idx, size_t col_idx, real_t val)
    {
        assert(row_idx < _row_dim);
        assert(col_idx < _col_dim);
        _entry_list.push_back(mat_entry{row_idx, col_idx, val});
    }
    void reserve(size_t entry_count) { _entry_list.reserve(entry_count); }
    size_t row_dim() const { return _row_dim; }
    size_t col_dim() const { return _col_dim; }
    size_t entry_count() const { return _entry_list.size(); }
    matrix create();
};

class sqr_mat_crtor : public matrix_crtor
{
public:
    sqr_mat_crtor(sqr_mat_crtor &&) = default;
    sqr_mat_crtor(size_t dim) : matrix_crtor(dim, dim) {}
    size_t dim() const { return _row_dim; }
    sqr_mat create();
};

inline real_t absolute(real_t v) { return v < 0 ? -v : v; }
}  // namespace markovgg
