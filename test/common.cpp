#include "common.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include "matrix.hpp"

using namespace markovgg;

void print_sep() { std::cout << "-------------------" << std::endl; }

spmatrix read_cord_format(std::istream& in)
{
    std::string line;
    size_t m = 0;
    size_t n = 0;
    size_t count = 0;
    while (!in.eof())
    {
        std::getline(in, line);
        if (line.size() > 0 && line[0] != '%')
        {
            std::istringstream sline(line);
            sline >> m;
            sline >> n;
            sline >> count;
            assert(!sline.fail());
            break;
        }
    }
    assert(m > 0);
    assert(n > 0);
    assert(count > 0);
    size_t i;
    size_t j;
    real_t v;
    spmatrix_creator crtor(m, n);
    while (!in.eof())
    {
        std::getline(in, line);
        if (line.size() > 0 && line[0] != '%')
        {
            std::istringstream sline(line);
            sline >> i;
            sline >> j;
            sline >> v;
            assert(!sline.fail());
            crtor.add_entry(i - 1, j - 1, v);
            count -= 1;
        }
    }
    assert(count == 0);
    return crtor.create(true);
}

matrix read_array_format(std::istream& in)
{
    std::string line;
    size_t m = 0;
    size_t n = 0;
    size_t count = 0;
    while (!in.eof())
    {
        std::getline(in, line);
        if (line.size() > 0 && line[0] != '%')
        {
            std::istringstream sline(line);
            sline >> m;
            sline >> n;
            assert(!sline.fail());
            break;
        }
    }
    assert(m > 0);
    assert(n > 0);
    real_t v;
    matrix A(m, n);
    while (!in.eof())
    {
        std::getline(in, line);
        if (line.size() > 0 && line[0] != '%')
        {
            std::istringstream sline(line);
            sline >> v;
            assert(!sline.fail());
            A(count % m, count / m) = v;
            count += 1;
        }
    }
    assert(count == m * n);
    return A;
}
