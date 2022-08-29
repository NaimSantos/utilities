#pragma once

#include <chrono>
#include <ctime>
#include <random>
#include <iostream>
#include <filesystem>

namespace ntime
{
	using milisegundos = std::chrono::milliseconds;
	using SteadyTimePoint = std::chrono::time_point<std::chrono::steady_clock>;

	class CustomTimer
	{
		private:
				SteadyTimePoint startpoint;
		public:
			CustomTimer(){
				startpoint = std::chrono::steady_clock::now();
			}
			~CustomTimer(){
				auto endpoint = std::chrono::steady_clock::now();
				auto start = std::chrono::time_point_cast<milisegundos>(startpoint).time_since_epoch().count();
				auto end = std::chrono::time_point_cast<milisegundos>(endpoint).time_since_epoch().count();
				auto total = end - start;

				std::cout << "\nElapsed time: " << total << " milliseconds" << std::endl;
			}
	};

	//Uma forma de capturar precisamente o tempo:
	SteadyTimePoint capture_time(){
		return std::chrono::steady_clock::now();;
	}

	//Retorna, em milissegundos, o tempo decorrido entre duas capturas de tempo:
	double get_elapsed_time(SteadyTimePoint ti, SteadyTimePoint tf){
		return static_cast<double>(std::chrono::duration_cast<milisegundos> (tf - ti).count());
	}

	//Exibe o horario atual em formato humanamente legível:
	void print_current_time(){
		std::time_t now = std::time(nullptr);
		std::cout << std::asctime(std::localtime(&now));
	}
}

//Obtains the names of all files in the current path, without checking for subpaths:
void GetFileNames(){
	for(auto& filepath: fs::directory_iterator(fs::current_path())){
		auto filename = filepath.path().filename();
		std::cout << filename << '\n';
	}
}



namespace nrandom
{
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
}
