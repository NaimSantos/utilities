#include "util_time.h"

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
