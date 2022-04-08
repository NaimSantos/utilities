#praga once

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