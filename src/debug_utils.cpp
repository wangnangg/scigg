#include "debug_utils.hpp"
#include "matrix.hpp"
namespace markovgg
{
void print_mat(const matrix& mat, uint8_t width)
{
    if (mat.row_dim() > 100 && mat.col_dim() > 100)
    {
        return;
    }
    std::cout << std::endl;
    for (size_t i = 0; i < mat.row_dim(); i++)
    {
        for (size_t j = 0; j < mat.col_dim(); j++)
        {
            std::cout << std::setw(width) << mat(i, j) << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void print_vec(const vec& v)
{
    if (v.dim() > 100)
    {
        return;
    }
    std::cout << "( ";
    for (size_t i = 0; i < v.dim(); i++)
    {
        std::cout << v[i] << " ";
    }
    std::cout << ")" << std::endl;
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
