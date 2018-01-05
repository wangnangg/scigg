#include "spmatrix.hpp"
#include <algorithm>
#include <utility>
#include "debug_utils.hpp"
namespace markovgg
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

real_t bin_search(spmat_vec_const_view view, size_t idx, real_t def_val)
{
    const size_t* begin = view.idx;
    const size_t* end = view.idx + view.nnz;
    const size_t* first = std::lower_bound(begin, end, idx);
    if (first != end && *first == idx)
    {
        return view.val[first - begin];
    }
    else
    {
        return def_val;
    }
}

real_t spmatrix_get_entry(spmatrix_const_view mat, size_t row, size_t col)
{
    if (mat.is_compressed_row())
    {
        auto view = mat[row];
        return bin_search(view, col, 0.0);
    }
    else
    {
        auto view = mat[col];
        return bin_search(view, row, 0.0);
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
