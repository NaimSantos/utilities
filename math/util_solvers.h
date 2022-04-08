#praga once

#include <vector>

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