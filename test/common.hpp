#pragma once
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include "matvec.hpp"
#include "spmatvec.hpp"
#include "type.hpp"

using namespace scigg;

void print_sep();

spmatrix read_cord_format(std::istream& in);
matrix read_array_format(std::istream& in);
