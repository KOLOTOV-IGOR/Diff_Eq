#include "Classes.h"

Solver::Solver() {
	Fill_M();
	complex<double> alpha(-0.5, -0.5);
	u = u0;
	result.assign(nt*nx, 0);
	for (int col = 0; col < nx; col++) {
		result[col] = u[col];
	}
	for (int row = 1; row < nt; row++) {
		Fill_Fu();
		Fill_F();
		k = Progonka(M_low + tau*alpha*Fu_low, M_diag + tau*alpha*Fu_diag, M_up + tau*alpha*Fu_up, F, k);
		u.erase(u.begin()); u.erase(u.end() - 1);
		u_new = u + tau*Real(k);
		for (int col = 1; col < nx - 1; col++) {
			result[row*nx + col] = u_new[col - 1];
		}
		u = u_new; u.insert(u.begin(), 0); u.insert(u.end(), 0);
		Fu_diag.clear(); Fu_up.clear(); Fu_low.clear();//Нужна отчистка, иначе проблемы 
		F.clear(); k.clear(); u_new.clear();
	}
	ofstream fout("result_new.txt");
	for (int row = 0; row < nt; row++) {
		for (int col = 0; col < (int)u0.size(); col++) {
			fout << result[row*(int)u0.size() + col] << " ";
		}
		fout << endl;
	}
	fout.close();
}

Solver::Solver(int t, double l, double pow_of_q, double time) {
	lambda = l; q = pow_of_q;
	Fill_M();
	complex<double> alpha(-0.5, -0.5);
	u = u0;
	result.assign(nt*nx, 0);
	for (int col = 0; col < nx; col++) {
		result[col] = u[col];
	}
	for (int row = 1; row < nt; row++) {
		Fill_Fu();
		Fill_F();
		k = Progonka(M_low + tau*alpha*Fu_low, M_diag + tau*alpha*Fu_diag, M_up + tau*alpha*Fu_up, F, k);
		u.erase(u.begin()); u.erase(u.end() - 1);
		u_new = u + tau*Real(k);
		for (int col = 1; col < nx - 1; col++) {
			result[row*nx + col] = u_new[col - 1];
		}
		u = u_new; u.insert(u.begin(), 0); u.insert(u.end(), 0);
		Fu_diag.clear(); Fu_up.clear(); Fu_low.clear();//Нужна отчистка, иначе проблемы 
		F.clear(); k.clear(); u_new.clear();
	}
	ofstream fout("result_new.txt");
	for (int row = 0; row < nt; row++) {
		for (int col = 0; col < (int)u0.size(); col++) {
			fout << result[row*(int)u0.size() + col] << " ";
		}
		fout << endl;
	}
}

Solver::~Solver() {
	M_diag.clear(); M_up.clear(); M_low.clear();
	Fu_diag.clear(); Fu_up.clear(); Fu_low.clear();//Нужна отчистка, иначе проблемы 
	F.clear(); k.clear(); u_new.clear();
}

int sign(double val) {
	if (val == 0) {
		return 0;
	}
	else if (val > 0) {
		return 1;
	}
	else
		return -1;
}

void Solver::Fill_M() {
	for (int i = 0; i < nx - 2; i++) {
		M_diag.push_back(lambda - (2 / pow(h, 2)));
		M_up.push_back(1 / pow(h, 2));
		M_low.push_back(1 / pow(h, 2));
	}
}

void Solver::Fill_F() {
	for (int i = 1; i < nx - 1; i++) {
		if (i == 1) {
			F.push_back((2 * u[i] - u[i + 1]) / pow(h, 2) - u[i] * pow(abs(u[i]), q - 1));
		}
		else if (i == nx - 2) {
			F.push_back((2 * u[i] - u[i - 1]) / pow(h, 2) - u[i] * pow(abs(u[i]), q - 1));
		}
		else {
			F.push_back(-(u[i + 1] - 2 * u[i] + u[i - 1]) / pow(h, 2) - u[i] * pow(abs(u[i]), q - 1));
		}
	}
}

void Solver::Fill_Fu() {
	for (int i = 1; i < nx - 1; i++) {
		Fu_diag.push_back(2 / pow(h, 2) - pow(abs(u[i]), q) - u[i] * (q - 1)*pow(abs(u[i]), q - 2)*sign(u[i]) - 1 / pow(h, 2));
		Fu_up.push_back(-1 / pow(h, 2));
		Fu_low.push_back(-1 / pow(h, 2));
	}
}

vec_complex operator+(vec_complex &a, vec_complex &b) {
	vec_complex temp;
	if (a.size() == b.size()) {
		for (int i = 0; i < (int)a.size(); i++) {
			temp.push_back(a[i] + b[i]);
		}
	}
	return temp;
}

vec operator+(vec &a, vec &b) {
	vec temp;
	if (a.size() == b.size()) {
		for (int i = 0; i < (int)a.size(); i++) {
			temp.push_back(a[i] + b[i]);
		}
	}
	return temp;
}

vec_complex operator*(complex<double> &a, vec_complex &b) {
	vector<complex<double>> temp;
	for (int i = 0; i < (int)b.size(); i++) {
		temp.push_back(a * b[i]);
	}
	return temp;
}

vec operator*(double &a, vec &b) {
	vector<double> temp;
	for (int i = 0; i < (int)b.size(); i++) {
		temp.push_back(a * b[i]);
	}
	return temp;
}

vec Real(vec_complex r) {
	vec temp;
	for (int i = 0; i < (int)r.size(); i++) {
		temp.push_back(real(r[i]));
	}
	return temp;
}

vec_complex Solver::Progonka(vec_complex &a, vec_complex &c, vec_complex &b, vec_complex &f, vec_complex &k) {
	int n = (int)c.size();
	vec_complex alpha; vec_complex beta;
	k.assign(n, 0);
	alpha.assign(n, 0); beta.assign(n, 0);
	alpha[0] = 0; beta[0] = 0;
	alpha[1] = -b[1] / c[1]; beta[1] = f[1] / c[1];
	for (int i = 1; i < n - 1; i++) {
		alpha[i + 1] = -b[i] / (a[i] * alpha[i] + c[i]);
		beta[i + 1] = (f[i] - a[i] * beta[i]) / (a[i] * alpha[i] + c[i]);
	}
	k[n - 1] = (f[n - 1] - a[n - 1] * beta[n - 1]) / (a[n - 1] * alpha[n - 1] + c[n - 1]);
	for (int j = n - 2; j >= 0; j--) {
		k[j] = alpha[j + 1] * k[j + 1] + beta[j + 1];
	}
	return k;
}
