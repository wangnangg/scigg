#include "matrix.hpp"
#include <algorithm>
#include <cassert>

namespace markovgg
{
bool cmp_vec_e(const sps_entry& x, const sps_entry& y) { return x.idx < y.idx; }

template <typename entry_type>
size_t count_sorted_unique(std::vector<entry_type>& list,
                           bool (*cmp)(const entry_type& x,
                                       const entry_type& y))
{
    size_t count = list.size();
    if (count == 0)
    {
        return count;
    }
    for (size_t i = 0; i < list.size() - 1; i++)
    {
        if (!cmp(list[i], list[i + 1]))
        {
            count -= 1;
        }
    }
    return count;
}

const sps_entry& bin_search(const sps_entry* first, const sps_entry* last,
                            const sps_entry& value,
                            bool (*cmp)(const sps_entry&, const sps_entry&))
{
    first = std::lower_bound(first, last, value, cmp);
    if (first != last && !cmp(value, *first))
    {
        return *first;
    }
    else
    {
        return value;
    }
}

bool cmp_mat_e(const mat_entry& e1, const mat_entry& e2)
{
    return e1.col_idx != e2.col_idx ? e1.col_idx < e2.col_idx
                                    : e1.row_idx < e2.row_idx;
}

bool cmp_sqr_mat_e(const mat_entry& e1, const mat_entry& e2)
{
    if (e1.col_idx < e2.col_idx)
    {
        return true;
    }
    if (e1.col_idx > e2.col_idx)
    {
        return false;
    }
    else
    {
        int64_t e1i = e1.row_idx == e1.col_idx ? -1 : e1.row_idx;
        int64_t e2i = e2.row_idx == e2.col_idx ? -1 : e2.row_idx;
        return e1i < e2i;
    }
}

matrix matrix_crtor::create()
{
    std::sort(_entry_list.begin(), _entry_list.end(), cmp_mat_e);
    size_t uniq_count = count_sorted_unique(_entry_list, cmp_mat_e);
    std::vector<size_t> idx;
    idx.reserve(_col_dim + 1);
    std::vector<sps_entry> val;
    val.reserve(uniq_count);
    idx.push_back(0);
    auto itor = _entry_list.begin();
    for (size_t i = 0; i < _col_dim; i++)
    {
        size_t counter = 0;
        while (itor != _entry_list.end() && itor->col_idx == i)
        {
            sps_entry entry{itor->row_idx, itor->val};
            itor += 1;
            while (itor != _entry_list.end() && itor->col_idx == i &&
                   itor->row_idx == entry.idx)
            {
                entry.val += itor->val;
                itor += 1;
            }
            if (entry.val != 0.0)
            {
                val.push_back(entry);
                counter += 1;
            }
        }
        idx.push_back(counter + idx[i]);
    }
    assert(_col_dim + 1 == idx.size());
    return matrix(_row_dim, std::move(idx), std::move(val));
}

real_t matrix::operator()(size_t row_idx, size_t col_idx) const
{
    const auto& col = this->col(col_idx);
    return bin_search(col.begin(), col.end(), sps_entry{row_idx, 0.0},
                      cmp_vec_e)
        .val;
}

sqr_mat sqr_mat_crtor::create()
{
    std::sort(_entry_list.begin(), _entry_list.end(), cmp_sqr_mat_e);
    size_t uniq_count = count_sorted_unique(_entry_list, cmp_sqr_mat_e);
    std::vector<size_t> idx;
    idx.reserve(dim() + 1);
    std::vector<sps_entry> val;
    val.reserve(uniq_count);
    idx.push_back(0);
    auto itor = _entry_list.begin();
    for (size_t i = 0; i < dim(); i++)
    {
        size_t counter = 0;
        if (itor == _entry_list.end() || itor->col_idx != i)
        {
            val.push_back({i, 0.0});
            counter += 1;
        }
        else
        {
            if (itor->row_idx == itor->col_idx)
            {
                sps_entry entry{itor->row_idx, itor->val};
                itor += 1;
                while (itor != _entry_list.end() && itor->col_idx == i &&
                       itor->row_idx == entry.idx)
                {
                    entry.val += itor->val;
                    itor += 1;
                }
                val.push_back(entry);
                counter += 1;
            }
            else
            {
                val.push_back({i, 0.0});
                counter += 1;
            }
            while (itor != _entry_list.end() && itor->col_idx == i)
            {
                sps_entry entry{itor->row_idx, itor->val};
                itor += 1;
                while (itor != _entry_list.end() && itor->col_idx == i &&
                       itor->row_idx == entry.idx)
                {
                    entry.val += itor->val;
                    itor += 1;
                }
                if (entry.val != 0.0)
                {
                    val.push_back(entry);
                    counter += 1;
                }
            }
        }
        idx.push_back(counter + idx[i]);
    }
    assert(dim() + 1 == idx.size());
    return sqr_mat(dim(), std::move(idx), std::move(val));
}

real_t sqr_mat::operator()(size_t row_idx, size_t col_idx) const
{
    if (row_idx == col_idx)
    {
        return this->diag(row_idx);
    }
    const auto& col = this->col(col_idx);
    return bin_search(col.begin() + 1, col.end(), sps_entry{row_idx, 0.0},
                      cmp_vec_e)
        .val;
}

bool operator==(const vec& x, const vec& y) { return x._val == y._val; }

bool operator==(const sps_entry& x, const sps_entry& y)
{
    return x.idx == y.idx && x.val == y.val;
}

}
