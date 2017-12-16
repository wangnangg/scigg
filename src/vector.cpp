#include "vector.hpp"

namespace markovgg
{
bool operator==(const vector &x, const vector &y) { return x._data == y._data; }
};
