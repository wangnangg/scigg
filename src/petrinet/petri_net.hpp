#pragma once
#include <algorithm>
#include <functional>
#include <vector>
#include "debug.hpp"
#include "type.hpp"

namespace scigg
{
template <typename Tdata, typename Marking, typename Token>
class petri_net
{
public:
    class context
    {
    private:
        const petri_net& pn;
        const Marking& mk;

    public:
        context(const petri_net& pn, const Marking& mk) : pn(pn), mk(mk) {}
        bool is_absorbing();
        Token token(size_t pid);
    };

    template <typename R>
    class marking_dep_var
    {
    public:
        using func_t = std::function<R(context ctxt)>;
        marking_dep_var() : _func(nullptr) {}
        marking_dep_var(func_t func) : _func(func) {}
        marking_dep_var(R value) : _val(value), _func(nullptr) {}
        marking_dep_var(R (*func_ptr)(context ctxt)) : _func(func_ptr) {}
        R operator()(context ctxt) const
        {
            if (_func)
            {
                return _func(ctxt);
            }
            else
            {
                return _val;
            }
        }
        bool is_val() { return !_func; }

    private:
        R _val;
        func_t _func;
    };
    using mk_uint = marking_dep_var<uint_t>;
    using mk_bool = marking_dep_var<bool>;

    struct arc
    {
        size_t pid;
        mk_uint multi;
        arc(size_t pid, mk_uint multi) : pid(pid), multi(multi) {}
    };

    struct transition
    {
        const petri_net* pn;
        size_t idx;
        uint_t prio;
        mk_bool enable;
        Tdata data;
        std::vector<arc> in_arc;
        std::vector<arc> out_arc;
        std::vector<arc> inh_arc;
        transition(const petri_net* pn, size_t idx, uint_t prio, mk_bool enable,
                   Tdata data)
            : pn(pn), idx(idx), prio(prio), enable(enable), data(data)
        {
        }
        transition& add_in_arc(size_t pid, mk_uint multi)
        {
            assert(pn->valid_pid(pid));
            in_arc.emplace_back(pid, multi);
            return *this;
        }
        transition& add_out_arc(size_t pid, mk_uint multi)
        {
            assert(pn->valid_pid(pid));
            out_arc.emplace_back(pid, multi);
            return *this;
        }
        transition& add_inh_arc(size_t pid, mk_uint multi)
        {
            assert(pn->valid_pid(pid));
            inh_arc.emplace_back(pid, multi);
            return *this;
        }
    };

private:
    static bool trans_comp(const transition& tr1, const transition& tr2)
    {
        return tr1.prio > tr2.prio;
    }
    std::vector<transition> _trans_list;
    std::vector<size_t> _trans_idx_map;
    size_t _place_count;
    mk_bool _g_enable;

public:  // observer
    std::vector<Marking> reached_markings(const Marking& m) const;
    const transition& operator[](size_t tid) const { return find_trans(tid); }
    std::vector<size_t> enabled_trans(const Marking& mk) const
    {
        auto etrans = std::vector<size_t>();
        bool found_enabled = false;
        uint_t enabled_prio;
        for (const auto& tr : _trans_list)
        {
            if (found_enabled && enabled_prio > tr.prio)
            {
                return etrans;
            }
            if (can_trans_enabled(tr, mk))
            {
                etrans.push_back(tr.idx);
                enabled_prio = tr.prio;
                found_enabled = true;
            }
        }
        return etrans;
    }
    bool is_trans_enabled(size_t t_ind, const Marking& mk) const
    {
        const auto& trans = find_trans(t_ind);
        for (const auto& tr : _trans_list)
        {
            if (tr.prio == trans.prio)
            {
                break;
            }
            if (can_trans_enabled(tr, mk))
            {
                return false;
            }
        }
        return can_trans_enabled(trans, mk);
    }

    Marking fire_trans(size_t t_ind, const Marking& mk) const
    {
        assert(is_trans_enabled(t_ind, mk));
        auto result_mk = mk;
        const auto& tr = find_trans(t_ind);
        context ct = {*this, mk};
        for (const auto& arc : tr.in_arc)
        {
            auto multi = arc.multi(ct);
            assert(result_mk[arc.pid] >= multi);
            result_mk[arc.pid] -= multi;
        }
        for (const auto& arc : tr.out_arc)
        {
            auto multi = arc.multi(ct);
            result_mk[arc.pid] += multi;
        }
        return result_mk;
    }

    Marking empty_marking() const { return Marking(_place_count, 0); }

public:  // modifier
    petri_net(size_t place_count) : _place_count(place_count), _g_enable(true)
    {
    }
    transition& add_trans(uint_t prio, mk_bool enable, Tdata data)
    {
        size_t new_idx = _trans_list.size();
        auto tr = transition(this, new_idx, prio, enable, std::move(data));
        auto hprio_tr = std::upper_bound(_trans_list.begin(), _trans_list.end(),
                                         tr, trans_comp);
        _trans_list.insert(hprio_tr, std::move(tr));
        _trans_idx_map.resize(_trans_list.size());
        for (size_t i = 0; i < _trans_list.size(); i++)
        {
            _trans_idx_map[_trans_list[i].idx] = i;
        }
        return find_trans(new_idx);
    }
    transition& add_in_arc(size_t tid, size_t pid, mk_uint multi);
    transition& add_out_arc(size_t tid, size_t pid, mk_uint multi);
    transition& add_inh_arc(size_t tid, size_t pid, mk_uint multi);
    void set_g_enable(mk_bool enable) { _g_enable = enable; }

private:
    transition& find_trans(size_t tid)
    {
        return _trans_list[_trans_idx_map[tid]];
    }
    const transition& find_trans(size_t tid) const
    {
        return _trans_list[_trans_idx_map[tid]];
    }

    bool can_trans_enabled(const transition& tr, const Marking& mk) const
    {
        context ct = {*this, mk};
        if (!_g_enable(ct))
        {
            return false;
        }
        for (const auto& arc : tr.in_arc)
        {
            if (mk[arc.pid] < arc.multi(ct))
            {
                return false;
            }
        }
        for (const auto& arc : tr.inh_arc)
        {
            if (mk[arc.pid] >= arc.multi(ct))
            {
                return false;
            }
        }
        if (!tr.enable(ct))
        {
            return false;
        }
        return true;
    }
    bool valid_pid(size_t pid) const { return pid < _place_count; }
};
}
