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
std::ostream& operator<<(std::ostream& os, const spmat_cs_entry& e);
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
}
