#pragma once
#include <cstdint>
#include <exception>
#include <string>
namespace markovgg
{
using std::size_t;
typedef double real_t;
typedef int int_t;
class msg_exception : public std::exception
{
public:
    msg_exception(std::string msg) : _msg(std::move(msg)) {}
    const char* what() const noexcept override { return _msg.c_str(); }

private:
    std::string _msg;
};
}
