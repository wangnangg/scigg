#pragma once
#include <cassert>
#include <vector>
#include "type.hpp"
namespace markovgg
{
enum spmatrix_format
{
    CPR_ROW,  // compressed row storage
    CPR_COL,  // compressed column storage
};
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

class spmatrix
{
    spmatrix_format _format;
    size_t _m;
    size_t _n;
    std::vector<size_t> _ptr;
    std::vector<spmat_cs_entry> _data;

public:
    spmatrix(spmatrix_format format, size_t m, size_t n,
             std::vector<size_t> ptr, std::vector<spmat_cs_entry> data)
        : _format(format),
          _m(m),
          _n(n),
          _ptr(std::move(ptr)),
          _data(std::move(data))
    {
        if (format == CPR_ROW)
        {
            assert(_ptr.size() - 1 == m);
        }
        else
        {
            assert(_ptr.size() - 1 == n);
        }
    }
    spmat_vec_const_view operator[](size_t i) const
    {
        return spmat_vec_const_view{&_data[_ptr[i]], &_data[_ptr[i + 1]]};
    }
    spmat_vec_mutable_view operator[](size_t i)
    {
        return spmat_vec_mutable_view{&_data[_ptr[i]], &_data[_ptr[i + 1]]};
    }
    real_t operator()(size_t row, size_t col) const;
    spmatrix_format format() const { return _format; }
    size_t m() const { return _m; }
    size_t n() const { return _n; }
    size_t ldim() const { return _ptr.size() - 1; }
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
    spmatrix_creator(size_t m, size_t n) : _m(m), _n(n), _data() {}
    void add_entry(size_t row, size_t col, real_t val)
    {
        _data.push_back(spmat_triplet_entry{row, col, val});
    }
    spmatrix create(spmatrix_format format);
};
}
