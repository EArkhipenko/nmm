#pragma once
#include <iostream>
#include<stdio.h>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#define EPS 1.e-30

using namespace std;

class Matrix
{
	typedef double real;
	typedef double real2;
	int * ia;
	real * di;
	real * al;

	vector<real> f;

	int n; // Размерность
	int m; // Кол-во элементов (размер al)
public:

	Matrix() {
		n = 0;
		m = 0;
		ia = nullptr;
	}


	int line_kol(int i, int m) {
		return i == 1 ? 1 : i % 2 + 2;
	}

	// Строим портрет матрицы по сетке
	void build(int kol, int koldof) {
		n = (kol + 1) * koldof / 2;
		di = new real[n];
		f.resize(n);
		ia = new int[n + 1];
		ia[0] = ia[1] = 0;
		for (int i = 2; i <= n; ++i) {
			ia[i] = ia[i - 1] + line_kol(i - 1, koldof);
		}
		
		m = ia[n];
		al = new real[m];
		for (int i = 0; i < n; ++i) {
			f[i] = di[i] = 0;
		}
		for (int i = 0; i < m; ++i) {
			al[i] = 0;
		}
	}

	void clean() {
		for (int i = 0; i < n; ++i) {
			f[i] = di[i] = 0;
		}
		for (int i = 0; i < m; ++i) {
			al[i] = 0;
		}
	}

	// Добавление локальной матрицы в глобальную
	void add(real koef, real * matr, int * ind, int k) {
		for (int i = 0; i < k; ++i) {
			di[ind[i]] += koef * matr[i * k + i]; // Добавили диагональ
			if (i == 0) { continue; }
			int ia0 = ia[ind[i]];
			int ia1 = ia[ind[i] + 1];
			int pj = ind[i] - (ia1 - ia0);

			for (int j = 0, s = ia0 + ind[j] - pj; s < ia1; ++s, ++j) {
				al[s] += koef * matr[i * k + j];
			}
		}
	}

	// Добавление локального вектора в глобальный
	void add_f(real koef, real * vec, int * ind, int k) {
		for (int i = 0; i < k; ++i) {
			f[ind[i]] += koef * vec[i];
		}
	}

	// Учет первых краевых условий
	void set_first(int n_dof, int dof, real rp) {
		di[dof] = 1.0; f[dof] = rp;
		if (dof == 0) {
			for (int i = 1; i < n_dof; ++i) {
				f[i] -= f[dof] * al[ia[i]];
				al[ia[i]] = 0.0;
			}
		}
		else if (dof == n - 1) {
			for (int i = n - 2, k = ia[dof + 1] - 1; i >= n - n_dof; --i, k--) {
				f[i] -= f[dof] * al[k];
				al[k] = 0.0;
			}
		}
		else {
			printf("Неверный номер узла для первого краевого условия");
		}
	}

	// Получить правую часть
	real & rp(int i) {
		return f[i];
	}

	// Получить диагональ
	real & diag(int i) {
		return di[i];
	}

	real element(int i, int j) {
		if (i == j) return di[i];
		if (j > i) std::swap(i, j); // Меняем местами элементы.
		int pi = i - ia[i + 1] + ia[i];
		if (j < pi) return 0;
		return al[ia[i] + j - pi];
	}

	void set_element(int i, int j, double val) {
		if (i == j) {
			di[i] = val;
			return;
		}
		if (j > i) std::swap(i, j); // Меняем местами элементы.
		int pi = i - ia[i + 1] + ia[i];
		if (j < pi) return;
	    al[ia[i] + j - pi] = val;
	}

	void print() {
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j) {
				printf("%f  ", element(i, j));
			}
			printf("| %f\n", f[i]);
		}
		//for (int i = 0; i < n; ++i) {
		//	printf("%f ", di[i]);
		//}printf("\n");
		//for (int i = 0; i < m; ++i) {
		//	printf("%f ", al[i]);
		//}
	}

	int size() {
		return n;
	}

	real sign(real a) {
		return a < 0 ? -1 : 1;
	}

	void LLt() {
		for (int i = 0; i < n; ++i) { // Цикл по строкам
			int ia0 = ia[i];
			int ia1 = ia[i + 1];
			real2 sd = 0;
			for (int k = ia0, j = i - (ia1 - ia0); k < ia1; ++k, ++j) { 
				real2 s = 0;

				int ja0 = ia[j];
				int ja1 = ia[j + 1];

				int pj = j - (ja1 - ja0);
				int pi = i - (ia1 - ia0);
				int ki = ia0;
				int kj = ja0;
				if (pj > pi) {
					ki += pj - pi;
				}
				else if (pi > pj) {
					kj += pi - pj;
				}

				for (; ki < k; ++ki, ++kj) {
					s += al[ki] * al[kj];
				}
				al[k] = (al[k] - s) / di[j];
				sd += al[k] * al[k];
			}
			// Вычисляем диагональный элемент
			real t = di[i] - sd;
			if (t < EPS) { 
				printf("di[%d] == %lf\n", i, t); 
				exit(1);
				return; 
			}
			di[i] = std::sqrt(t);
		}
	}

	void pryamoi_hod(const vector<real> & f, vector<real> & x) {
		for (int i = 0; i < n; ++i) {
			real2 s = 0;
			int ia0 = ia[i];
			int ia1 = ia[i + 1];
			for (int k = ia0, j = i - (ia1 - ia0); k < ia1; ++k, ++j) {
				s += al[k] * x[j];
			}
			if (di[i] < EPS) {
				printf("di[%d] = %f < EPS\n", i, di[i]);
				exit(1);
			}
			x[i] = (f[i] - s) / di[i];
		}
	}

	void obratniy_hod(vector<real> & f, vector<real> & x) {
		for (int j = n-1; j >= 0; j--) {
			if (di[j] < EPS) {
				printf("di[%d] = %f < EPS\n", j, di[j]);
				exit(1);
			}
			x[j] = f[j] / di[j];
			int k = ia[j + 1] - ia[j];
			for (int l = ia[j], i = j - k; l < ia[j + 1]; ++l, ++i) {
				f[i] -= x[j] * al[l];
			}
		}
	}

	void solve(vector<real> & x) {
		vector<real> z(size());
		LLt();
		pryamoi_hod(f, z);
		obratniy_hod(z, x);
	}
};





