#pragma once

#include <random>
#include <chrono>

int random_int(int lower_index, int upper_index);
template <typename T = double> T random_floatpoint(T lower_index, T upper_index);

