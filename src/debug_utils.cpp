#include "debug_utils.hpp"
#include "matrix.hpp"
namespace markovgg
{
std::ostream& operator<<(std::ostream& os, const matrix_const_view& mat)
{
    if (mat.m() > 100 && mat.n() > 100)
    {
        return os << "matrix (" << mat.m() << ", " << mat.n() << ")";
    }
    for (size_t i = 0; i < mat.m(); i++)
    {
        for (size_t j = 0; j < mat.n(); j++)
        {
            os << std::setw(6) << mat(i, j) << " ";
        }
        os << std::endl;
    }
    return os;
}

std::ostream& operator<<(std::ostream& os, const vector_const_view& v)
{
    if (v.dim() > 100)
    {
        return os << "vector (" << v.dim() << ")";
    }
    os << "[ ";
    for (size_t i = 0; i < v.dim(); i++)
    {
        os << v[i] << " ";
    }
    return os << "]" << std::endl;
}
timed_scope::timed_scope(const char* name)
{
    std::cout << name << " timer started." << std::endl;
    _start_time = std::chrono::high_resolution_clock::now();
    _name = name;
}

timed_scope::~timed_scope()
{
    std::chrono::duration<real_t> time_span =
        std::chrono::duration_cast<std::chrono::duration<real_t>>(
            std::chrono::high_resolution_clock::now() - _start_time);
    std::cout << _name << " used " << time_span.count() << " seconds"
              << std::endl;
}
}
