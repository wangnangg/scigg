#pragma once
#include <cassert>
#include <vector>
#include "type.hpp"
namespace markovgg
{
struct spmat_cs_entry
{
    size_t idx;
    real_t val;
};

struct spmat_vec_const_view  // make for loop happy
{
    const spmat_cs_entry* const entry_start;
    const spmat_cs_entry* const entry_end;

    const spmat_cs_entry* begin() const { return entry_start; }
    const spmat_cs_entry* end() const { return entry_end; }
};

struct spmat_vec_mutable_view  // make for loop happy, do not ever change idx
                               // field
{
    spmat_cs_entry* const entry_start;
    spmat_cs_entry* const entry_end;

    spmat_cs_entry* begin() const { return entry_start; }
    spmat_cs_entry* end() const { return entry_end; }
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
    const size_t* _ptr;
    spmat_cs_entry* _data;
    friend class spmatrix_const_view;
    friend class spmatrix;

public:
    spmatrix_mutable_view(size_t m, size_t n, bool is_crs, const size_t* ptr,
                          spmat_cs_entry* data)
        : spmatrix_base(m, n, is_crs), _ptr(ptr), _data(data)
    {
    }
    spmatrix_mutable_view(spmatrix& mat);
    spmat_vec_mutable_view operator[](size_t i) const
    {
        return spmat_vec_mutable_view{&_data[_ptr[i]], &_data[_ptr[i + 1]]};
    }
    real_t operator()(size_t row, size_t col) const;
    spmatrix_mutable_view transpose() const
    {
        return spmatrix_mutable_view(n(), m(), !is_compressed_row(), _ptr,
                                     _data);
    }
};

class spmatrix_const_view : public spmatrix_base
{
    const size_t* _ptr;
    const spmat_cs_entry* _data;
    friend class spmatrix;

public:
    spmatrix_const_view(size_t m, size_t n, bool is_crs, const size_t* ptr,
                        const spmat_cs_entry* data)
        : spmatrix_base(m, n, is_crs), _ptr(ptr), _data(data)
    {
    }
    spmatrix_const_view(const spmatrix& mat);
    spmatrix_const_view(spmatrix_mutable_view mat)
        : spmatrix_base(mat.m(), mat.n(), mat.is_compressed_row()),
          _ptr(mat._ptr),
          _data(mat._data)
    {
    }

    spmat_vec_const_view operator[](size_t i) const
    {
        return spmat_vec_const_view{&_data[_ptr[i]], &_data[_ptr[i + 1]]};
    }
    real_t operator()(size_t row, size_t col) const;
    spmatrix_const_view transpose() const
    {
        return spmatrix_const_view(n(), m(), !is_compressed_row(), _ptr, _data);
    }
};

class spmatrix : public spmatrix_base
{
    std::vector<size_t> _ptr;
    std::vector<spmat_cs_entry> _data;
    friend class spmatrix_mutable_view;
    friend class spmatrix_const_view;

public:
    spmatrix(size_t m, size_t n, bool is_row_compressed,
             std::vector<size_t> ptr, std::vector<spmat_cs_entry> data)
        : spmatrix_base(m, n, is_row_compressed),
          _ptr(std::move(ptr)),
          _data(std::move(data))
    {
        if (is_row_compressed)
        {
            assert(_ptr.size() - 1 == m);
        }
        else
        {
            assert(_ptr.size() - 1 == n);
        }
    }
    spmatrix_const_view transpose() const
    {
        return spmatrix_const_view(n(), m(), !is_compressed_row(), &_ptr[0],
                                   &_data[0]);
    }
    spmatrix_mutable_view transpose()
    {
        return spmatrix_mutable_view(n(), m(), !is_compressed_row(), &_ptr[0],
                                     &_data[0]);
    }
    spmat_vec_mutable_view operator[](size_t i)
    {
        return spmat_vec_mutable_view{&_data[_ptr[i]], &_data[_ptr[i + 1]]};
    }
    spmat_vec_const_view operator[](size_t i) const
    {
        return spmat_vec_const_view{&_data[_ptr[i]], &_data[_ptr[i + 1]]};
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
