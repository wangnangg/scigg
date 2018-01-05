#pragma once
#include <cstdint>
#include <exception>
#include <string>
namespace markovgg
{
using std::size_t;
typedef double real_t;
typedef int int_t;
typedef unsigned int uint_t;
struct msg_exception : public std::exception
{
    std::string msg;
    msg_exception() = default;
    msg_exception(std::string msg) : msg(std::move(msg)) {}
    const char* what() const noexcept override { return msg.c_str(); }
};
}
