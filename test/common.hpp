#pragma once
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include "matrix.hpp"
#include "spmatrix.hpp"
#include "type.hpp"

using namespace markovgg;

void print_sep();

spmatrix read_cord_format(std::istream& in);
matrix read_array_format(std::istream& in);
