#include "matrix_oper.hpp"
#include <cassert>

namespace markovgg
{
bool near_eq(const vec &x, const vec &y, double prec)
{
    for (size_t i = 0; i < x.dim(); i++)
    {
        if (!near_zero(x[i] - y[i], prec))
        {
            return false;
        }
    }
    return true;
}
bool near_zero(const vec &x, double prec)
{
    for (size_t i = 0; i < x.dim(); i++)
    {
        if (!near_zero(x[i], prec))
        {
            return false;
        }
    }
    return true;
}
void fill(vec &x, real_t val)
{
    for (size_t i = 0; i < x.dim(); i++)
    {
        x[i] = val;
    }
}

real_t norm1(const vec &vec)
{
    real_t norm = 0;
    for (auto e : vec)
    {
        norm += absolute(e);
    }
    return norm;
}

void set_norm1(vec &vec, real_t val)
{
    real_t norm = norm1(vec);
    real_t factor = val / norm;
    for (auto &e : vec)
    {
        e *= factor;
    }
}

real_t vec_sum(const vec &v)
{
    real_t vsum = 0;
    for (auto e : v)
    {
        vsum += e;
    }
    return vsum;
}

void set_vec_sum(vec &v, real_t val)
{
    real_t vsum = vec_sum(v);
    real_t factor = val / vsum;
    for (auto &e : v)
    {
        e *= factor;
    }
}

void acc_sub(vec &x, const vec &v, const vec &w)
{
    for (size_t i = 0; i < x.dim(); i++)
    {
        x[i] += v[i] - w[i];
    }
}
void sub(vec &x, const vec &v, const vec &w)
{
    for (size_t i = 0; i < x.dim(); i++)
    {
        x[i] = v[i] - w[i];
    }
}
void acc_add(vec &x, const vec &v, const vec &w)
{
    for (size_t i = 0; i < x.dim(); i++)
    {
        x[i] += v[i] + w[i];
    }
}
void add(vec &x, const vec &v, const vec &w)
{
    for (size_t i = 0; i < x.dim(); i++)
    {
        x[i] = v[i] + w[i];
    }
}

void self_add(vec &x, const vec &v)
{
    for (size_t i = 0; i < x.dim(); i++)
    {
        x[i] += v[i];
    }
}
void acc_mul(vec &x, const vec &v, real_t a)
{
    for (size_t i = 0; i < x.dim(); i++)
    {
        x[i] += v[i] * a;
    }
}
void mul(vec &x, const vec &v, real_t a)
{
    for (size_t i = 0; i < x.dim(); i++)
    {
        x[i] = v[i] * a;
    }
}
void acc_mul(vec &x, const vec &v, const matrix &A)
{
    for (size_t i = 0; i < A.col_dim(); i++)
    {
        for (const auto &entry : A.col(i))
        {
            x[i] += entry.val * v[entry.idx];
        }
    }
}

void mul(vec &x, const vec &v, const matrix &A)
{
    fill(x, 0.0);
    acc_mul(x, v, A);
}

void scale(matrix &A, real_t a)
{
    for (size_t i = 0; i < A.col_dim(); i++)
    {
        for (auto &e : A.col(i))
        {
            e.val *= a;
        }
    }
}

void scale(vec &v, real_t a)
{
    for (auto &e : v)
    {
        e *= a;
    }
}

const size_t *search_new_idx(const std::vector<size_t> &idx_list, size_t idx)
{
    auto found_iter = std::lower_bound(idx_list.begin(), idx_list.end(), idx);
    if (found_iter != idx_list.end() && *found_iter == idx)
    {
        return &*found_iter;
    }
    else
    {
        return nullptr;
    }
}

sqr_mat part_mat(const matrix &mat, const std::vector<size_t> &sorted_valid_idx)
{
    sqr_mat_crtor crtor(sorted_valid_idx.size());
    for (size_t col = 0; col < sorted_valid_idx.size(); col++)
    {
        for (auto &e : mat.col(sorted_valid_idx[col]))
        {
            auto idx_ptr = search_new_idx(sorted_valid_idx, e.idx);
            if (idx_ptr != nullptr)
            {
                crtor.add_entry(idx_ptr - &*sorted_valid_idx.begin(), col,
                                e.val);
            }
        }
    }
    return crtor.create();
}

matrix part_mat(const matrix &mat, const std::vector<size_t> &sorted_row_idx,
                const std::vector<size_t> &sorted_col_idx)
{
    size_t row_dim = sorted_row_idx.size();
    size_t col_dim = sorted_col_idx.size();
    matrix_crtor crtor(row_dim, col_dim);
    for (size_t col = 0; col < sorted_col_idx.size(); col++)
    {
        for (auto &e : mat.col(sorted_col_idx[col]))
        {
            auto idx_ptr = search_new_idx(sorted_row_idx, e.idx);
            if (idx_ptr != nullptr)
            {
                crtor.add_entry(idx_ptr - &*sorted_row_idx.begin(), col, e.val);
            }
        }
    }
    return crtor.create();
}

vec part_vec(const vec &v, const std::vector<size_t> &sorted_valid_idx)
{
    vec p(sorted_valid_idx.size(), 0.0);
    for (size_t i = 0; i < p.dim(); i++)
    {
        p[i] = v[sorted_valid_idx[i]];
    }
    return p;
}

sqr_mat transpose(const sqr_mat &mat)
{
    sqr_mat_crtor crtor(mat.dim());
    for (size_t i = 0; i < mat.dim(); i++)
    {
        for (auto &e : mat.col(i))
        {
            crtor.add_entry(i, e.idx, e.val);
        }
    }
    return crtor.create();
}
}
