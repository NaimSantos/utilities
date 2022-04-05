#pragma once
#include <chrono>
#include <ctime>
#include <iostream>

// Objeto de registro automatico de tempo através de escopo:
class CustomTimer;

// Captura precisa do tempo atual:
std::chrono::time_point<std::chrono::steady_clock> capture_time();
// Retorna, em milissegundos, o tempo decorrido entre duas capturas de tempo:
double get_elapsed_time(std::chrono::time_point<std::chrono::steady_clock> ti, std::chrono::time_point<std::chrono::steady_clock> tf);
// Horario atual em formato humanamente legivel:
void print_current_time();