#include <fstream>
#include "common.hpp"
#include "debug_utils.hpp"
#include "gtest/gtest.h"
TEST(test_read, read_cord_format)
{
    {
        std::ifstream file("./test/matrix/example/example.mtx");
        spmatrix mat = read_cord_format(file);
        print(mat);
        file.close();
    }
    {
        std::ifstream file("./test/matrix/example/example_b.mtx");
        matrix mat = read_array_format(file);
        print(mat);
        file.close();
    }
}
