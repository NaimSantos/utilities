#pragma once

#include <vector>
#include <iostream>
#include <functional>
#include <cmath>

/* Parâmetros:
xi = valor inicial de x
xf = valor final de x
y0 = valor inicial da função
fxy = a função que é avaliada
*/

double solveEuler(double xi, double xf, double y0, double h, std::function<double (double, double)> fxy){
	auto n = std::round((xf - xi) / h);
	if (n <=0 )
		return NAN;

	for (int i = 1; i <= n; i++){
		y0 = y0 + fxy(xi, y0)*h;
		xi = xi + h;
	}
	return y0;
}

double solveEuler_Mel(double xi, double xf, double y0, double h, std::function<double (double, double)> fxy){
	auto n = std::round((xf - xi) / h);
	if (n <=0 )
		return NAN;

	double k1{0}, k2{0};
	double y = y0;
	
	for (int i = 1; i <= n; i++){
		k1 = fxy(xi, y);
		k2 = fxy(xi + h, y + k1*h);

		y = y + (k1 + k2)*0.5*h;
		xi = xi + h;
	}
	return y;
}

double solveEuler_Mod(double xi, double xf, double y0, double h, std::function<double (double, double)> fxy){
	auto n = std::round((xf - xi) / h);
	if (n <=0 )
		return NAN;

	double k1{0}, k2{0};
	double y = y0;

	for (size_t i {1}; i <= n; i++){
		k1 = fxy(xi, y);
		k2 = fxy(xi + 0.5*h, y + 0.5*k1*h);

		y = y + k2*h;
		xi = xi + h;
	}
	return y;
}

double solveHuen_C(double xi, double xf, double y0, double h, std::function<double (double, double)> fxy){
	auto n = std::round((xf - xi) / h);
	if (n <=0 )
		return NAN;

	double k1{0}, k2{0};
	double y = y0;

	for (size_t i {1}; i <= n; i++){
		k1 = fxy(xi, y);
		k2 = fxy(xi + h, y + k1*h);

		y = y + (0.5*k1 + 0.5*k2)*h;
		xi = xi + h;
	}
	return y;
}

double solveRalston(double xi, double xf, double y0, double h, std::function<double (double, double)> fxy){
	auto n = std::round((xf - xi) / h);
	if (n <=0 )
		return NAN;

	double k1{0}, k2{0};
	double y = y0;

	for (size_t i {1}; i <= n; i++){
		k1 = fxy(xi, y);
		k2 = fxy(xi + 0.75*h, y + 0.75*k1*h);

		y = y + ((1.0/3.0)*k1 + (2.0/3.0)*k2)*h;
		xi = xi + h;
	}
	return y;
}

double solveRK3(double xi, double xf, double y0, double h, std::function<double (double, double)> fxy){
	auto n = std::round((xf - xi) / h);
	if (n <=0 )
		return NAN;

	double k1{0}, k2{0}, k3{0};
	double y = y0;

	for (size_t i {1}; i <=n ; i++){
		k1 = fxy(xi, y);
		k2 = fxy(xi + 0.5*h, y + 0.5*k1*h);
		k3 = fxy(xi + h, y - k1*h + 2*k2*h);

		y = y + (1.0/6.0)*(k1 + 4*k2 + k3)*h;
		xi = xi + h;
	}
	return y;
}

double solveRK4(double xi, double xf, double y0, double h, std::function<double (double, double)> fxy){
	auto n = std::round((xf - xi) / h);
	if (n <=0 )
		return NAN;

	double k1{0}, k2{0}, k3{0}, k4{0};
	double y = y0;

	for (size_t i {1}; i <= n ; i++){
		k1 = fxy(xi, y);
		k2 = fxy(xi + 0.5*h, y + 0.5*k1*h);
		k3 = fxy(xi + 0.5*h, y + 0.5*k2*h);
		k4 = fxy(xi + h, y + k3*h);

		y = y + (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4)*h;
		xi = xi + h;
	}
	return y;
}

double solveRK5(double xi, double xf, double y0, double h, std::function<double (double, double)> fxy){
	auto n = std::round((xf - xi) / h);
	if (n <=0 )
		return NAN;

	double k1{0}, k2{0}, k3{0}, k4{0}, k5{0}, k6{0};
	double y = y0;

	for (size_t i {1}; i <=n ; i++){
		k1 = fxy(xi, y);
		k2 = fxy(xi + 0.25*h, y + 0.25*k1*h);
		k3 = fxy(xi + 0.25*h, y + 0.125*k1*h + 0.125*k2*h);
		k4 = fxy(xi + 0.5*h, y - 0.5*k2*h + k3*h);
		k5 = fxy(xi + 0.75*h, y + (3.0/16.0)*k1*h + (9.0/16.0)*k4*h);
		k6 = fxy(xi + h, y - (3.0/7.0)*k1*h + (2.0/7.0)*k2*h + (12.0/7.0)*k3*h - (12.0/7.0)*k4*h + (8.0/7.0)*k5*h);

		y = y + (1.0/90.0)*(7*k1 + 32*k2 + 12*k4 + 32*k5 + 7*k6)*h;
		xi = xi + h;
	}
	return y;
}



void solveAdams(const double h, double xi, const double xf, std::vector<double>& Y, std::function<double (double, double)> fxy){

	solveRK4(xi, xf-h, Y[0], h, fxy); //Runge kutta determina os valores não fornecidos
	const auto n = std::round((xf - xi) / h);

	double milne_p {0.0}, milne_c {0.0};
	auto x = xi;
	for (size_t i = 4; i <= n; ++i){

		auto f4 = fxy(x+i*h, Y[i]);
		auto f3 = fxy(x+(i-1)*h, Y[i-1]);
		auto f2 = fxy(x+(i-2)*h, Y[i-2]);
		auto f1 = fxy(x+(i-3)*h, Y[i-3]);

		milne_p = Y[i] + ((55.0/24.0)*f4 - (59.0/24.0)*f3 + (37.0/24.0)*f2 - (9.0/24.0)*f1)*h;
		milne_c = Y[i] + ((9.0/24.0)*fxy(x+i*h, milne_p) + (19.0/24.0)*fxy(x+(i)*h, Y[i]) - (5.0/24.0)*fxy(x+(i-1)*h, Y[i-1]) + (1.0/24.0)*fxy(x+(i-2)*h, Y[i-2]))*h;

		Y.push_back(milne_c);
	}
}

// Técnicas para obtenção de raizes de equações:
double solveBisect(double a0, double b0, std::function<double (double)> ftarget){
	if (ftarget(a0) * ftarget(b0) >= 0){
		return NAN;
	}
	double c0 {a0};
	double eps {0.00001};
	unsigned int i{1};
	while ((b0-a0) >= eps || i<100){
		c0 = (a0 + b0) / 2;
		if (ftarget(c0) == 0)
			break;
		else
			(ftarget(c0)*ftarget(a0) < 0) ? (b0 = c0) : (a0 = c0);
		i++;
	}
	if(i>=100){
		return NAN;
	}
	return c0;
}

//A derivada é avaliada por diferencas finitas centrada. Erro: O(h^4)
double f_diff(double x, double h, std::function<double (double)> f){
	double res = (-f(x+2*h) + 8*f(x+h) - 8*f(x-h) + f(x-2*h)) / (12*h);
	return res;
}

double solveNewton(double x, std::function <double (double)> fnwtn){
	auto h = fnwtn(x) / f_diff(x, 0.001, fnwtn);
	unsigned int i {0};
	while (std::fabs(h) >= 0.00001 && i<100){
		h = fnwtn(x) / f_diff(x, 0.001, fnwtn);
		x = x - h;
		i++;
	}
	return x;
}

//Diferenças finitas avançada para primeira derivada.
double fdiff_a_1st(double x, double h, std::function<double (double)> f){
	return (f(x+h) - f(x)) / h;
}

double fdiff_a_1st_v2(double x, double h, std::function<double (double)> f){
	return (-f(x+2*h) + 4*f(x+h) - 3*f(x)) / (2*h);
}

//Diferenças finitas avançada para segunda derivada.
double fdiff_a_2nd(double x, double h, std::function<double (double)> f){
	return (f(x+2*h) - 2*f(x+h) + f(x)) / (h*h);
}
double fdiff_a_2nd_v2(double x, double h, std::function<double (double)> f){
	return (-f(x+3*h) - 4*f(x+2*h) - 5*f(x+h) + 2*f(x)) / (h*h);
}

//Diferenças finitas recuada para primeira derivada.
double fdiff_r_1st(double x, double h, std::function<double (double)> f){
	return (f(x) - f(x-h)) / h;
}
double fdiff_r_1st_v2(double x, double h, std::function<double (double)> f){
	return (3*f(x) - 4*f(x-h) + f(x-2*h)) / (2*h);
}

//Diferenças finitas recuada para segunda derivada:
double fdiff_r_2nd(double x, double h, std::function<double (double)> f){
	return (f(x) - 2*f(x-h) + f(x-2*h)) / (h*h);
}
double fdiff_r_2nd_v2(double x, double h, std::function<double (double)> f){
	return (2*f(x) - 5*f(x-h) + 4*f(x-2*h) + f(x-3*h)) / (h*h);
}

//Diferenças finitas centrada  para primeira derivada:
double fdiff_c_1st(double x, double h, std::function<double (double)> f){
	return (f(x+h) - f(x-h)) / (2*h);
}
double fdiff_c_1st_v2(double x, double h, std::function<double (double)> f){
	return (-f(x+2*h) + 8*f(x+h) - 8*f(x-h) + f(x-2*h)) / (12*h);
}

//Diferenças finitas centrada  para segunda derivada:
double fdiff_c_2nd(double x, double h, std::function<double (double)> f){
	return (f(x+h) - 2*f(x) + f(x-h)) / (h*h);
}
double fdiff_c_2nd_v2(double x, double h, std::function<double (double)> f){
	return (-f(x+2*h) + 16*f(x+h) - 30*f(x) + 16*f(x-h) - f(x-2*h)) / (12*(h*h));
}

// Técnicas de integração numérica:
double int_trapz(double a, double b, double h, std::function<double (double)> f){
	double res {0.0};
	const auto n = static_cast<int>( std::floor((std::fabs(b - a)) / h));
	for (int i = 0; i < n - 1; i++){
		res += f( a + i*h);
	}
	res += (f(a) + f(b) ) / 2;
	res *= h;

	return res;
}
double int_simpson(double a, double b, double h, std::function<double (double)> f){	
	const auto n = static_cast<int>( std::floor((std::fabs(b - a)) / h));
	//n - 1 pontos internos
	double sum_odds {0.0};
	for (int i = 1; i < n; i += 2){
		sum_odds += f(a + i*h);
	}
	double sum_evens {0.0};
	for (int i = 2; i < n; i += 2){
		sum_evens += f(a + i*h);
	}
	return (f(a) + f(b) + 2*sum_evens + 4*sum_odds) * (h/3);
}


/*Resolve um sistema utilizando o TDMA
a = diagonal inferior (std::vector<double>)
b = diagonal principal (std::vector<double>)
c = diagonal superior (std::vector<double>)
d = vetor de termos independentes

c e d são modificados. d conterá a solução
*/
void solve_by_tdma(const std::vector<double>& a, const std::vector<double>& b, std::vector<double>& c, std::vector<double>& d){
	auto n = static_cast<int>(d.size()-1);

	c[0] /= b[0];
	d[0] /= b[0];

	for (int i = 1; i < n; i++){
		c[i] = (c[i] ) / (b[i] - a[i]*c[i-1]);
		d[i] = (d[i] - a[i]*d[i-1]) / (b[i] - a[i]*c[i-1]);
	}

	d[n] = (d[n] - a[n]*d[n-1]) / (b[n] - a[n]*c[n-1]);

	for (int i = n; i-- > 0;){
		d[i] = d[i] - (c[i]*d[i+1]);
	}
}

/* Versão alternativa, que não altera c e d, mas requer um vetor adicional f, para armazenar a solução
*/
void solve_by_tdma(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c,  const std::vector<double>& d, std::vector<double>& f){
	int N = d.size();

	std::vector<double> c_start(N, 0.0);
	std::vector<double> d_start(N, 0.0);

	c_start[0] = c[0] / b[0];
	d_start[0] = d[0] / b[0];

	for (int i=1; i<N; i++){
		double m = 1.0 / (b[i] - a[i] * c_start[i-1]);
		c_start[i] = c[i] * m;
		d_start[i] = (d[i] - a[i] * d_start[i-1]) * m;
	}

	for (int i=N-1; i-- > 0; ) {
		f[i] = d_start[i] - c_start[i] * d[i+1];
	}
}

//Uma estimativa para o valor de PI:
constexpr double NPI = 4 * std::atan(1);
//Conversão de graus para radianos:
double to_rad(double value){
	return value * NPI / 180;
}

//Retorna se um arranjo bidimensional é uma matriz diagonalmente dominante:
template <typename T>
bool is_diagonal_dom(T** M, const int nrows, const int ncol){
	if (nrows != ncol)
		return false;

	double sum {0.0};
	for (int i = 0; i < nrows; i++){
		double diag = std::fabs(M[i][i]);
		for (int j = 0; j < ncol; j++){
			if (i == j)
				continue;
			sum += std::fabs(M[i][j]);
		}
		if (sum > diag)
			return false;
	}
	return true;
}

std::vector<double> linspace(const int n, const double xf, const double xi){
	auto h = (xf - xi) / (n-1);
	std::vector<double> Vec (n, 0.0);
	for (int i = 0; i < n; i++){
		Vec[i] = xi + i*h;
	}
	return Vec;
}

struct LinReg
{
	double a {0.0};
	double b {0.0};
	double r {0.0};
	double eps {0.0};
};

//Funções auxiliares requeridas pela regressão:
template <typename T>
double average(T* A, const unsigned int n){
	double sum {0.0};
	for (int i = 0; i < n ; i++)
		sum += A[i];
	return sum / n;
}

template <typename T>
double desvio(T* A, const unsigned int n, const T aver){
	double sum {0.0};
	for (int i = 0; i < n; i++){
		sum += std::pow(A[i] - aver, 2) / (n - 1);
	}
	return std::sqrt(sum);
}

template <typename T1, typename T2>
double covar(T1* A1, T2* A2, const unsigned int n){
	auto x_med = average<T1>(A1, n);
	auto y_med = average<T2>(A2, n);
	double sum {0.0};

	for (int i = 0; i < n; i++){
		sum += ((A1[i] - x_med)*(A2[i] - y_med)) / (n - 1);
	}
	return sum;
}

/*
Retorna um struct com 4 componentes, coeficientes da regressão:
	a: inclinação da reta
	b: interseção da reta
	r: coeficiente de correlação de Pearson
	eps: erro padrão em y
O parâmetro de template T1 é o tipo dos dados de 'x'; T2 é o tipo dos dados de 'y'.
*/
template<typename T1, typename T2>
LinReg RegressaoLinear(T1* A1, T2* A2, const unsigned int n){

	LinReg res;

	auto covarxy = covar<T1, T2>(A1, A2, n);
	auto xmed = average<T1>(A1, n);
	auto ymed = average<T2>(A2, n);
	auto desvx = desvio<T1>(A1, n, xmed);
	auto desvy = desvio<T2>(A2, n, ymed);

	auto a = (covarxy / (desvx*desvx));
	res.a = a;

	auto b = ymed - (a * xmed);
	res.b = b;

	auto r = covarxy /( desvx * desvy);
	res.r = r;

	double sum = 0.0;
	for (int i = 0; i < n; i++){
		sum += std::pow((A2[i] - ((a*A1[i]) + b)), 2) / (n - 2);
	}
	auto eps = std::sqrt(sum);
	res.eps = eps;

	return res;
}
