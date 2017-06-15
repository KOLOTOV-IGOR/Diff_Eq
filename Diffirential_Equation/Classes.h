#pragma once
#ifndef CLASSES_H
#define CLASSES_H
#include<iostream>
#include<vector>
#include<iterator>
#include<math.h>
#include<complex>
#include<fstream>
#include<time.h> 
using namespace std;
typedef vector<complex<double>> vec_complex;
typedef vector<double> vec;

class Grid {
protected:
	double rb, lb;//rb - right bound, rb - right bound
	double T;//upper time bound
	double ulb, urb;//lb - left bound, rb - right bound
	double h, tau; //steps
	int nx, nt;
	vec u0;//Начальные данные
public:
	Grid() {
		lb = 0, rb = 3.141592, T = 1, nt = 50;
		ulb = -8; urb = 4;
		u0.clear();
		ifstream f("file_for_u0.txt");
		while (!f.eof()) {
			double v;
			f >> v;
			u0.push_back(v);
		}
		f.close();
		u0.erase(u0.end() - 1);//Без этого будет дублироваться последнее значение.
		nx = (int)u0.size();
		get_h(); get_tau();
	}
	Grid(int t, double time) {
		lb = 0, rb = 3.141592, T = time, nt = t;
		urb = 0; ulb = 0;
		u0.clear();
		ifstream f("file_for_u0.txt");
		while (!f.eof()) {
			double v;
			f >> v;
			u0.push_back(v);
		}
		f.close();
		u0.erase(u0.end() - 1);//Без этого будет дублироваться последнее значение.
		nx = (int)u0.size();
		get_h(); get_tau();
	}
	~Grid() { }
	void get_h() {
		h = (rb - lb) / nx;
	}
	void get_tau() {
		tau = T / nt;
	}
	void get_bound_condition(double left, double right) {
		ulb = left; urb = right;
	}
};

class Solver : public Grid {
private:
	double lambda, q;
	vec u;//Решение на предыдущем шаге.
	vec u_new;//Решение на следующем шаге.
	vec result;
	/////////Матрицы для Розенброка
	vec_complex F;
	vec_complex k;
	vec_complex M_diag, M_up, M_low;
	vec_complex Fu_diag;
	vec_complex Fu_up;
	vec_complex Fu_low;
	//////////////////////////////////////////
public:
	Solver();
	Solver(int t, double l, double pow_of_q, double time);
	~Solver();
	void Fill_M();
	void Fill_Fu();
	void Fill_F();
	vec_complex Progonka(vec_complex &a, vec_complex &c, vec_complex &b, vec_complex &f, vec_complex &k);
	friend vec Real(vec_complex r);
	friend vec_complex operator+(vec_complex &a, vec_complex &b);
	friend vec_complex operator*(complex<double> &a, vec_complex &b);
	friend vec operator+(vec &a, vec &b);
	friend vec operator*(double &a, vec &b);
};

#endif CLASSES_H