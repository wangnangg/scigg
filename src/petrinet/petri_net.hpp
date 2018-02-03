#pragma once
#include <algorithm>
#include <functional>
#include <vector>
#include "debug.hpp"
#include "type.hpp"

namespace scigg
{
template <typename Tdata, typename Marking>
class petri_net_tmpl
{
public:
    using pid_t = size_t;
    using tid_t = size_t;
    using marking = Marking;
    using token = typename Marking::token;
    using tdata = Tdata;
    class context
    {
    private:
        const petri_net_tmpl& pn;
        const marking& mk;

    public:
        context(const petri_net_tmpl& pn, const marking& mk) : pn(pn), mk(mk) {}
        token ntoken(pid_t pid) { return mk[pid]; }
    };

    template <typename R>
    class marking_dep_var
    {
    public:
        using query_func = std::function<R(context ctxt)>;
        marking_dep_var() : _func(nullptr) {}
        template <typename QueryFunc,
                  std::enable_if_t<!std::is_convertible<QueryFunc, R>::value,
                                   int> = 0>
        marking_dep_var(QueryFunc func) : _func(func)
        {
        }
        marking_dep_var(R value) : _val(value), _func(nullptr) {}
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
        query_func _func;
    };
    using mk_uint = marking_dep_var<uint_t>;
    using mk_bool = marking_dep_var<bool>;
    using mk_real = marking_dep_var<real_t>;
    using mk_tdata = marking_dep_var<tdata>;

    class arc
    {
        pid_t _pid;
        mk_uint _multi;

    public:
        arc(pid_t pid, mk_uint multi) : _pid(pid), _multi(multi) {}
        pid_t pid() const { return _pid; }
        uint multi(context ct) const { return _multi(ct); };
    };

    class transition
    {
        const petri_net_tmpl* _pn;
        tid_t _idx;
        uint_t _prio;
        mk_bool _enable;
        mk_tdata _data;
        std::vector<arc> _in_arc;
        std::vector<arc> _out_arc;
        std::vector<arc> _inh_arc;

    public:
        transition(const petri_net_tmpl* pn, tid_t idx, uint_t prio,
                   mk_bool enable, mk_tdata data)
            : _pn(pn), _idx(idx), _prio(prio), _enable(enable), _data(data)
        {
        }
        transition& add_in_arc(pid_t pid, mk_uint multi = 1)
        {
            assert(_pn->valid_pid(pid));
            _in_arc.emplace_back(pid, multi);
            return *this;
        }
        transition& add_out_arc(pid_t pid, mk_uint multi = 1)
        {
            assert(_pn->valid_pid(pid));
            _out_arc.emplace_back(pid, multi);
            return *this;
        }
        transition& add_inh_arc(pid_t pid, mk_uint multi = 1)
        {
            assert(_pn->valid_pid(pid));
            _inh_arc.emplace_back(pid, multi);
            return *this;
        }
        tid_t idx() const { return _idx; }
        uint_t prio() const { return _prio; }
        const std::vector<arc>& in_arc() const { return _in_arc; };
        const std::vector<arc>& out_arc() const { return _out_arc; }
        const std::vector<arc>& inh_arc() const { return _inh_arc; };
        bool enabled(context ct) const { return _enable(ct); }
        tdata data(context ct) const { return _data(ct); }
    };

private:
    static bool trans_comp(const transition& tr1, const transition& tr2)
    {
        return tr1.prio() > tr2.prio();
    }
    std::vector<transition> _trans_list;
    std::vector<tid_t> _trans_idx_map;
    pid_t _place_count;
    mk_bool _g_enable;

public:  // observer
    std::vector<marking> reached_markings(const marking& m) const;
    const transition& operator[](tid_t tid) const { return find_trans(tid); }
    std::vector<tid_t> enabled_trans(const marking& mk) const
    {
        auto etrans = std::vector<tid_t>();
        bool found_enabled = false;
        uint_t enabled_prio;
        for (const auto& tr : _trans_list)
        {
            if (found_enabled && enabled_prio > tr.prio())
            {
                return etrans;
            }
            if (can_trans_enable(tr, mk))
            {
                etrans.push_back(tr.idx());
                enabled_prio = tr.prio();
                found_enabled = true;
            }
        }
        return etrans;
    }
    bool is_trans_enabled(tid_t t_ind, const marking& mk) const
    {
        const auto& trans = find_trans(t_ind);
        for (const auto& tr : _trans_list)
        {
            if (tr.prio() == trans.prio())
            {
                break;
            }
            if (can_trans_enable(tr, mk))
            {
                return false;
            }
        }
        return can_trans_enable(trans, mk);
    }

    marking fire_trans(tid_t t_ind, const marking& mk) const
    {
        assert(is_trans_enabled(t_ind, mk));
        auto result_mk = mk;
        const auto& tr = find_trans(t_ind);
        context ct = {*this, mk};
        for (const auto& arc : tr.in_arc())
        {
            auto multi = arc.multi(ct);
            assert(result_mk[arc.pid()] >= multi);
            result_mk[arc.pid()] -= multi;
        }
        for (const auto& arc : tr.out_arc())
        {
            auto multi = arc.multi(ct);
            result_mk[arc.pid()] += multi;
        }
        return result_mk;
    }

    tdata trans_data(tid_t t_ind, const marking& mk) const
    {
        return find_trans(t_ind).data({*this, mk});
    }

    marking empty_marking() const { return marking(_place_count, 0); }

public:  // modifier
    petri_net_tmpl(pid_t place_count)
        : _place_count(place_count), _g_enable(true)
    {
    }
    transition& add_trans(uint_t prio, mk_bool enable, mk_tdata data)
    {
        tid_t new_idx = _trans_list.size();
        auto tr = transition(this, new_idx, prio, enable, data);
        auto hprio_tr = std::upper_bound(_trans_list.begin(), _trans_list.end(),
                                         tr, trans_comp);
        _trans_list.insert(hprio_tr, std::move(tr));
        _trans_idx_map.resize(_trans_list.size());
        for (tid_t i = 0; i < _trans_list.size(); i++)
        {
            _trans_idx_map[_trans_list[i].idx()] = i;
        }
        return find_trans(new_idx);
    }
    transition& add_in_arc(tid_t tid, pid_t pid, mk_uint multi = 1)
    {
        return this->operator[](tid).add_in_arc(pid, multi);
    }
    transition& add_out_arc(tid_t tid, pid_t pid, mk_uint multi = 1)
    {
        return this->operator[](tid).add_out_arc(pid, multi);
    }
    transition& add_inh_arc(tid_t tid, pid_t pid, mk_uint multi = 1)
    {
        return this->operator[](tid).add_inh_arc(pid, multi);
    }
    void set_g_enable(mk_bool enable) { _g_enable = enable; }

private:
    transition& find_trans(tid_t tid)
    {
        return _trans_list[_trans_idx_map[tid]];
    }
    const transition& find_trans(tid_t tid) const
    {
        return _trans_list[_trans_idx_map[tid]];
    }

    bool can_trans_enable(const transition& tr, const marking& mk) const
    {
        context ct = {*this, mk};
        if (!_g_enable(ct))
        {
            return false;
        }
        for (const auto& arc : tr.in_arc())
        {
            if (mk[arc.pid()] < arc.multi(ct))
            {
                return false;
            }
        }
        for (const auto& arc : tr.inh_arc())
        {
            if (mk[arc.pid()] >= arc.multi(ct))
            {
                return false;
            }
        }
        if (!tr.enabled(ct))
        {
            return false;
        }
        return true;
    }
    bool valid_pid(pid_t pid) const { return pid < _place_count; }
};

template <typename Token>
class marking_tmpl
{
public:
    using token = Token;

private:
    std::vector<token> _token_list;

public:
    marking_tmpl(size_t ntoken, token val) : _token_list(ntoken, val) {}
    token operator[](size_t i) const { return _token_list[i]; }
    token& operator[](size_t i) { return _token_list[i]; }
    size_t size() const { return _token_list.size(); }
};

template <typename Marking>
struct marking_comp_tmpl
{
    bool operator()(const Marking& m1, const Marking& m2) const
    {
        for (size_t i = 0; i < m1.size(); i++)
        {
            if (m1[i] < m2[i])
            {
                return true;
            }
            else if (m1[i] > m2[i])
            {
                return false;
            }
        }
        return false;
    }
};

using default_marking = marking_tmpl<uint_t>;
using default_marking_comp = marking_comp_tmpl<default_marking>;

template <typename Tdata>
using default_petri_net = petri_net_tmpl<Tdata, default_marking>;
}
