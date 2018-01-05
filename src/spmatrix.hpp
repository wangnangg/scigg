#pragma once
#include <cassert>
#include <vector>
#include "type.hpp"

namespace markovgg
{
struct spmat_vec_const_view
{
    size_t nnz;
    const size_t* idx;
    const real_t* val;
};

struct spmat_vec_mutable_view
{
    size_t nnz;
    size_t* idx;
    real_t* val;
};

class spmatrix_base
{
protected:
    bool _is_crs;
    size_t _m;
    size_t _n;

public:
    spmatrix_base(size_t m, size_t n, bool is_compressed_row)
        : _is_crs(is_compressed_row), _m(m), _n(n)
    {
    }
    bool is_compressed_row() const { return _is_crs; }
    size_t m() const { return _m; }
    size_t n() const { return _n; }
    size_t ldim() const
    {
        if (is_compressed_row())
        {
            return m();
        }
        else
        {
            return n();
        }
    }
};

class spmatrix;

class spmatrix_mutable_view : public spmatrix_base
{
    size_t* _ptr;
    size_t* _idx;
    real_t* _val;
    friend class spmatrix_const_view;
    friend class spmatrix;

public:
    spmatrix_mutable_view(size_t m, size_t n, bool is_crs, size_t* ptr,
                          size_t* idx, real_t* val)
        : spmatrix_base(m, n, is_crs), _ptr(ptr), _idx(idx), _val(val)
    {
    }
    spmatrix_mutable_view(spmatrix& mat);
    spmat_vec_mutable_view operator[](size_t i) const
    {
        return spmat_vec_mutable_view{_ptr[i + 1] - _ptr[i], &_idx[_ptr[i]],
                                      &_val[_ptr[i]]};
    }
    real_t operator()(size_t row, size_t col) const;
    spmatrix_mutable_view transpose() const
    {
        return spmatrix_mutable_view(n(), m(), !is_compressed_row(), _ptr, _idx,
                                     _val);
    }
};

class spmatrix_const_view : public spmatrix_base
{
    const size_t* _ptr;
    const size_t* _idx;
    const real_t* _val;
    friend class spmatrix;

public:
    spmatrix_const_view(size_t m, size_t n, bool is_crs, const size_t* ptr,
                        const size_t* idx, const real_t* val)
        : spmatrix_base(m, n, is_crs), _ptr(ptr), _idx(idx), _val(val)
    {
    }
    spmatrix_const_view(const spmatrix& mat);
    spmatrix_const_view(spmatrix_mutable_view mat)
        : spmatrix_base(mat.m(), mat.n(), mat.is_compressed_row()),
          _ptr(mat._ptr),
          _idx(mat._idx),
          _val(mat._val)
    {
    }

    spmat_vec_const_view operator[](size_t i) const
    {
        return spmat_vec_const_view{_ptr[i + 1] - _ptr[i], &_idx[_ptr[i]],
                                    &_val[_ptr[i]]};
    }
    real_t operator()(size_t row, size_t col) const;
    spmatrix_const_view transpose() const
    {
        return spmatrix_const_view(n(), m(), !is_compressed_row(), _ptr, _idx,
                                   _val);
    }
};

class spmatrix : public spmatrix_base
{
    std::vector<size_t> _ptr;
    std::vector<size_t> _idx;
    std::vector<real_t> _val;
    friend class spmatrix_mutable_view;
    friend class spmatrix_const_view;

public:
    spmatrix(size_t m, size_t n, bool is_row_compressed,
             std::vector<size_t> ptr, std::vector<size_t> idx,
             std::vector<real_t> val)
        : spmatrix_base(m, n, is_row_compressed),
          _ptr(std::move(ptr)),
          _idx(std::move(idx)),
          _val(std::move(val))
    {
        if (is_row_compressed)
        {
            assert(_ptr.size() - 1 == m);
        }
        else
        {
            assert(_ptr.size() - 1 == n);
        }
        assert(_idx.size() == _val.size());
    }
    spmatrix_const_view transpose() const
    {
        return spmatrix_const_view(n(), m(), !is_compressed_row(), &_ptr[0],
                                   &_idx[0], &_val[0]);
    }
    spmatrix_mutable_view transpose()
    {
        return spmatrix_mutable_view(n(), m(), !is_compressed_row(), &_ptr[0],
                                     &_idx[0], &_val[0]);
    }
    spmat_vec_mutable_view operator[](size_t i)
    {
        return spmat_vec_mutable_view{_ptr[i + 1] - _ptr[i], &_idx[_ptr[i]],
                                      &_val[_ptr[i]]};
    }
    spmat_vec_const_view operator[](size_t i) const
    {
        return spmat_vec_const_view{_ptr[i + 1] - _ptr[i], &_idx[_ptr[i]],
                                    &_val[_ptr[i]]};
    }
    real_t operator()(size_t row, size_t col) const;
};

struct spmat_triplet_entry
{
    size_t row;
    size_t col;
    real_t val;
};
class spmatrix_creator
{
    std::vector<spmat_triplet_entry> _data;
    size_t _m;
    size_t _n;

public:
    spmatrix_creator(size_t m, size_t n) : _data(), _m(m), _n(n) {}
    void add_entry(size_t row, size_t col, real_t val)
    {
        _data.push_back(spmat_triplet_entry{row, col, val});
    }
    spmatrix create(bool is_row_compressed);
};
}
