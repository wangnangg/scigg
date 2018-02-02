#pragma once
#include "autodiff.hpp"
#include "petri_net.hpp"

namespace scigg
{
template <typename Marking>
using diff1_spn = petri_net_tmpl<const diff1_expr*, Marking>;

using default_diff1_spn = diff1_spn<default_marking>;
}
