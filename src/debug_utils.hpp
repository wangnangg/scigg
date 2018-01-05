#pragma once
#include <cassert>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <vector>
#include "matrix.hpp"
#include "spmatrix.hpp"
#include "vector.hpp"
namespace markovgg
{
std::ostream& operator<<(std::ostream& os, const matrix& mat);
std::ostream& operator<<(std::ostream& os, const matrix_mutable_view& mat);
std::ostream& operator<<(std::ostream& os, const matrix_const_view& mat);
std::ostream& operator<<(std::ostream& os, const vector& v);
std::ostream& operator<<(std::ostream& os, const vector_const_view& v);
std::ostream& operator<<(std::ostream& os, const vector_mutable_view& v);
std::ostream& operator<<(std::ostream& os, const spmatrix& mat);
std::ostream& operator<<(std::ostream& os, const spmatrix_const_view& mat);
std::ostream& operator<<(std::ostream& os, const spmatrix_mutable_view& mat);
std::ostream& operator<<(std::ostream& os, const spmat_triplet_entry& e);
template <class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v)
{
    os << "(";
    for (typename std::vector<T>::const_iterator ii = v.begin(); ii != v.end();
         ++ii)
    {
        os << " " << *ii;
    }
    os << ")";
    return os;
}

template <typename M>
void print(const M& v)
{
    std::cout << v << std::endl;
}

class timed_scope
{
    std::chrono::high_resolution_clock::time_point _start_time;
    const char* _name;

public:
    timed_scope(const char* name);
    ~timed_scope();
};

class timer
{
    std::chrono::duration<real_t> _duration =
        std::chrono::duration<real_t>::zero();
    std::chrono::high_resolution_clock::time_point _start_time;
    const char* _name;

public:
    timer(const char* name) { _name = name; }
    void start() { _start_time = std::chrono::high_resolution_clock::now(); }
    void stop()
    {
        _duration += std::chrono::duration_cast<std::chrono::duration<real_t>>(
            std::chrono::high_resolution_clock::now() - _start_time);
        _start_time = std::chrono::high_resolution_clock::now();
    }
    ~timer()
    {
        std::cout << _name << " used " << _duration.count() << " seconds"
                  << std::endl;
    }
};
}
