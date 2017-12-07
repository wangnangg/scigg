#pragma once
#include <cassert>
#include <chrono>
#include <iomanip>
#include <iostream>

namespace markovgg
{
template <typename val_t>
void print(const val_t& s)
{
    std::cout << s << std::endl;
}

class matrix;
void print_mat(const matrix& mat, uint8_t width = 6);

class vec;
void print_vec(const vec& v);

class timed_scope
{
    std::chrono::high_resolution_clock::time_point _start_time;
    const char* _name;

public:
    timed_scope(const char* name);
    ~timed_scope();
};
}
