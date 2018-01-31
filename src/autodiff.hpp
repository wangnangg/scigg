#pragma once
#include <memory>
#include "matvec.hpp"
namespace scigg
{
class diff1_expr
{
public:
    virtual real_t eval(vector_const_view x,
                        vector_mutable_view grad) const = 0;
};

class diff1_var : public diff1_expr
{
    size_t _var_idx;

public:
    diff1_var(size_t var_idx) : _var_idx(var_idx) {}
    real_t eval(vector_const_view x, vector_mutable_view grad) const override;
};

class diff1_con : public diff1_expr
{
    real_t _val;

public:
    diff1_con(real_t val) : _val(val) {}
    real_t eval(vector_const_view x, vector_mutable_view grad) const override;
};

class diff1_sum : public diff1_expr
{
    std::vector<const diff1_expr*> _children;

public:
    diff1_sum() = default;
    template <typename... Args>
    diff1_sum(Args... exprs) : _children{exprs...}
    {
    }
    void add_child(const diff1_expr* child) { _children.push_back(child); }
    real_t eval(vector_const_view x, vector_mutable_view grad) const override;
};

// grad(f)*g*h + f*grad(g)*h + f*g*grad(h)
class diff1_prod : public diff1_expr
{
    std::vector<const diff1_expr*> _children;

public:
    diff1_prod() = default;
    template <typename... Args>
    diff1_prod(Args... exprs) : _children{exprs...}
    {
    }
    void add_child(const diff1_expr* child) { _children.push_back(child); }
    real_t eval(vector_const_view x, vector_mutable_view grad) const override;
};

class diff1_expr_store
{
    std::vector<std::unique_ptr<diff1_expr>> _expr_list;

public:
    template <typename T, typename... Args>
    T* custom(Args&&... args)
    {
        _expr_list.emplace_back(new T(std::forward<Args>(args)...));
        return static_cast<T*>(_expr_list.back().get());
    }
    diff1_var* var(size_t var_idx) { return custom<diff1_var>(var_idx); }
    diff1_con* con(real_t val) { return custom<diff1_con>(val); }
    template <typename... Args>
    diff1_sum* sum(Args... args)
    {
        return custom<diff1_sum>(args...);
    }
    template <typename... Args>
    diff1_prod* prod(Args... args)
    {
        return custom<diff1_prod>(args...);
    }
};
};
