#include "util_random.h"

int random_int(int lower_index, int upper_index)
{
	std::default_random_engine rd{static_cast<long unsigned int>(std::chrono::high_resolution_clock::now().time_since_epoch().count()) };
	std::mt19937 gen(rd());
	std::uniform_int_distribution<int> distrib(lower_index, upper_index);

	return distrib(gen);
}

template <typename T = double>
T random_floatpoint(T lower_index, T upper_index)
{
	std::default_random_engine rd{static_cast<long unsigned int>(std::chrono::high_resolution_clock::now().time_since_epoch().count()) };
	std::mt19937 gen(rd());
	std::uniform_real_distribution<T> distrib(lower_index, upper_index);

	return distrib(gen);
}