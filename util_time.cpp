#include "util_time.h"

using SteadyClock = std::chrono::steady_clock
using milisegundos = std::chrono::milliseconds;


class CustomTimer
{
	public:
		CustomTimer(){
			startpoint = std::chrono::steady_clock::now();
		}
		~CustomTimer(){
			TimeDestructor();
		}

		void TimeDestructor(){
			auto endpoint = std::chrono::steady_clock::now();

			auto start = std::chrono::time_point_cast<milisegundos>(startpoint).time_since_epoch().count();
			auto end = std::chrono::time_point_cast<milisegundos>(endpoint).time_since_epoch().count();

			auto total = end - start;
			double s = total*0.001;
			
			std::cout << "\nTempo decorrido:" << s << " segundos (" << total << " milisegundos)" << std::endl;
		}
	
		private:
			std::chrono::time_point<SteadyClock> startpoint;
};

//Uma forma de capturar precisamente o tempo:
std::chrono::time_point<SteadyClock> capture_time(){
	auto res = std::chrono::steady_clock::now();
	return res;
}

//Retorna, em milissegundos, o tempo decorrido entre duas capturas de tempo:
double get_elapsed_time(std::chrono::time_point<SteadyClock> ti, std::chrono::time_point<SteadyClock> tf){
	double elapsed_time_ms = std::chrono::duration_cast<milisegundos> (tf - ti).count();
	return elapsed_time_ms;
}

void print_current_time(){
	std::time_t now = std::time(nullptr);
	std::cout << std::asctime(std::localtime(&now));
}