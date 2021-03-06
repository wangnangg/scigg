#include "spmatrix.hpp"
#include <algorithm>
#include <utility>
#include "debug.hpp"
namespace scigg
{
spmatrix_mutable_view::spmatrix_mutable_view(spmatrix& mat)
    : spmatrix_base(mat.m(), mat.n(), mat.is_compressed_row()),
      _ptr(&mat._ptr[0]),
      _idx(&mat._idx[0]),
      _val(&mat._val[0])
{
}

spmatrix_const_view::spmatrix_const_view(const spmatrix& mat)
    : spmatrix_base(mat.m(), mat.n(), mat.is_compressed_row()),
      _ptr(&mat._ptr[0]),
      _idx(&mat._idx[0]),
      _val(&mat._val[0])
{
}

spmatrix::spmatrix(spmatrix_const_view view, bool is_row_compressed)
    : spmatrix_base(view.m(), view.n(), is_row_compressed),
      _ptr(is_row_compressed ? view.m() + 1 : view.n() + 1),
      _idx(view.nnz()),
      _val(view.nnz())
{
    if (view.is_compressed_row() == is_row_compressed)
    {
        std::copy(view._ptr, view._ptr + _ptr.size(), &_ptr[0]);
        std::copy(view._idx, view._idx + _idx.size(), &_idx[0]);
        std::copy(view._val, view._val + _val.size(), &_val[0]);
    }
    else
    {
        std::vector<size_t> counter(ldim(), 0);
        for (size_t i = 0; i < view.ldim(); i++)
        {
            auto vec = view[i];
            for (size_t j = 0; j < vec.nnz; j++)
            {
                size_t idx = vec.idx[j];
                counter[idx] += 1;
            }
        }
        size_t pos = 0;
        _ptr[0] = pos;
        for (size_t i = 0; i < counter.size(); i++)
        {
            pos += counter[i];
            _ptr[i + 1] = pos;
        }
        for (size_t i = 0; i < view.ldim(); i++)
        {
            auto vec = view[i];
            for (size_t j = 0; j < vec.nnz; j++)
            {
                size_t idx = vec.idx[j];
                real_t val = vec.val[j];
                size_t offset = _ptr[idx + 1] - counter[idx];
                counter[idx] -= 1;
                _val[offset] = val;
                _idx[offset] = i;
            }
        }
    }
}
real_t spmatrix_get_entry(spmatrix_const_view mat, size_t row, size_t col)
{
    if (mat.is_compressed_row())
    {
        auto view = mat[row];
        auto entry_ptr = search_spvec_entry(view, col);
        if (entry_ptr)
        {
            return *entry_ptr;
        }
        else
        {
            return 0.0;
        }
    }
    else
    {
        auto view = mat[col];
        auto entry_ptr = search_spvec_entry(view, row);
        if (entry_ptr)
        {
            return *entry_ptr;
        }
        else
        {
            return 0.0;
        }
    }
}

real_t spmatrix::operator()(size_t row, size_t col) const
{
    return spmatrix_get_entry(*this, row, col);
}
real_t spmatrix_mutable_view::operator()(size_t row, size_t col) const
{
    return spmatrix_get_entry(*this, row, col);
}
real_t spmatrix_const_view::operator()(size_t row, size_t col) const
{
    return spmatrix_get_entry(*this, row, col);
}

static bool cmp_CPR_ROW(const spmat_triplet_entry& e1,
                        const spmat_triplet_entry& e2)
{
    if (e1.row < e2.row)
    {
        return true;
    }
    else if (e1.row > e2.row)
    {
        return false;
    }
    else if (e1.col < e2.col)
    {
        return true;
    }
    else
    {
        return false;
    }
}

static bool cmp_CPR_COL(const spmat_triplet_entry& e1,
                        const spmat_triplet_entry& e2)
{
    if (e1.col < e2.col)
    {
        return true;
    }
    else if (e1.col > e2.col)
    {
        return false;
    }
    else if (e1.row < e2.row)
    {
        return true;
    }
    else
    {
        return false;
    }
}

spmatrix spmatrix_creator::create(bool is_row_compressed)
{
    if (is_row_compressed)
    {
        std::sort(_data.begin(), _data.end(), cmp_CPR_ROW);
        std::vector<size_t> ptr_list;
        ptr_list.reserve(_m + 1);
        ptr_list.push_back(0);
        std::vector<size_t> idx_list;
        idx_list.reserve(_data.size());
        std::vector<real_t> val_list;
        val_list.reserve(_data.size());
        auto itor = _data.begin();
        for (size_t i = 0; i < _m; i++)
        {
            size_t counter = 0;
            while (itor != _data.end() && itor->row == i)
            {
                size_t idx = itor->col;
                real_t val = itor->val;
                itor += 1;
                while (itor != _data.end() && itor->row == i &&
                       itor->col == idx)
                {
                    val += itor->val;
                    itor += 1;
                }
                if (val != 0.0)
                {
                    idx_list.push_back(idx);
                    val_list.push_back(val);
                    counter += 1;
                }
            }
            ptr_list.push_back(counter + ptr_list[i]);
        }
        assert(_m + 1 == ptr_list.size());
        return spmatrix(_m, _n, is_row_compressed, std::move(ptr_list),
                        std::move(idx_list), std::move(val_list));
    }
    else
    {
        std::sort(_data.begin(), _data.end(), cmp_CPR_COL);
        std::vector<size_t> ptr_list;
        ptr_list.reserve(_n + 1);
        ptr_list.push_back(0);
        std::vector<size_t> idx_list;
        idx_list.reserve(_data.size());
        std::vector<real_t> val_list;
        val_list.reserve(_data.size());
        auto itor = _data.begin();
        for (size_t i = 0; i < _n; i++)
        {
            size_t counter = 0;
            while (itor != _data.end() && itor->col == i)
            {
                size_t idx = itor->row;
                real_t val = itor->val;
                itor += 1;
                while (itor != _data.end() && itor->col == i &&
                       itor->row == idx)
                {
                    val += itor->val;
                    itor += 1;
                }
                if (val != 0.0)
                {
                    idx_list.push_back(idx);
                    val_list.push_back(val);
                    counter += 1;
                }
            }
            ptr_list.push_back(counter + ptr_list[i]);
        }
        assert(_n + 1 == ptr_list.size());
        return spmatrix(_m, _n, is_row_compressed, std::move(ptr_list),
                        std::move(idx_list), std::move(val_list));
    }
}
}
