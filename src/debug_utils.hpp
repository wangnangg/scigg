#pragma once
#include <cassert>
#include <chrono>
#include <iomanip>
#include <iostream>
#include "matrix.hpp"
#include "vector.hpp"
namespace markovgg
{
std::ostream& operator<<(std::ostream& os, const matrix& mat);
std::ostream& operator<<(std::ostream& os, const vector& v);

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
