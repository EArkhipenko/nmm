#pragma once
#include<stdio.h>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <set>
#include <complex>
#include <functional>
#define EPS 1.e-30
#define DEL( X ) if(X) delete [] X;

class Matrix
{
	typedef double block;
	bool is_bad(const block & val) {
		return val < 1.e-13;
	}

	typedef double real;
	typedef double real2;
	// Матрица
	int * ig;
	int * jg;
	block * di;
	block * al;
	block * au;
	block * pr;

	// Разложение
	block * ggl;
	block * ggu;
	block * gdi;


	int n; // Размерность
	int m; // Кол-во элементов (размер al)
	int maxiter;
public:

	Matrix() {
		n = 0;
		m = 0;
		ig = nullptr;
		jg = nullptr;
		al = nullptr;
		au = nullptr;
		di = nullptr;
	}


	void setMaxiter(int mit) {
		maxiter = mit;
	}

	// Строим портрет матрицы по сетке
	void build(int kol, int kol_dofs, int kol_local_dofs,  std::function<void(int, int *)> getind) {
		std::vector<std::set<int>> vlist(kol_dofs);   //vlist массив списков смежных вершин для кэ
		int nelem = kol;  //количество элементов в сетве
		int npsi = kol_local_dofs;  //количество базисных функций
		int * ind = new int[npsi];  //глобальные индексы базисных функций для элемента
		for (int elem = 0; elem < nelem; ++elem) {
			getind(elem, ind);
			for (int i = 0; i < npsi; ++i) {
				int gi = ind[i];
				for (int j = i + 1; j < npsi; ++j) {
					int gj = ind[j];
					if (gi > gj) std::swap(gi, gj);
					vlist[gj].insert(gi);
				}
			}
		}

		n = kol_dofs;
		ig = new int[n + 1];
		ig[0] = ig[1] = 0;
		for (int i = 1; i < n; ++i) {
			ig[i + 1] = ig[i] + vlist[i].size();
		}

		m = ig[n];
		jg = new int[m];

		for (int i = 0, k = 0; i < n; ++i) {
			for (auto & j : vlist[i]) {
				jg[k] = j; 
				k++;
			}
		}

		di = new block[n];
		pr = new block[n];
		al = new block[m];
		au = new block[m];
		
		clear();
	}

	void clear() {
		for (int i = 0; i < n; ++i) {
			di[i] = pr[i] = 0;
		}

		for (int i = 0; i < m; ++i) {
			al[i] = au[i] = 0;
		}
	}

	// Добавление локальной матрицы в глобальную
	void add(real koef, block * matr, int * ind, int k) {
		for (int i = 0; i < k; ++i) {
			di[ind[i]] += koef * matr[i * k + i]; // Добавили диагональ
			int ia0 = ig[ind[i]];
			int ia1 = ig[ind[i] + 1];
			for (int j = 0, s = ia0; s < ia1; ++s) {
				if (jg[s] == ind[j]) {
					al[s] += koef * matr[i * k + j];
					au[s] += koef * matr[j * k + i];
					++j;
				}
			}
		}
	}

	// Добавление локального вектора в глобальный
	void add_f(real koef, block * vec, int * ind, int k) {
		for (int i = 0; i < k; ++i) {
			pr[ind[i]] += koef * vec[i];
		}
	}

	// Учет первых краевых условий
	void set_first(int i, block val) {                   //i строчка в которой задаём первое краевое  val  значение которое ставим в правую часть
		di[i] = 1.0; 
		pr[i] = val;
		for (int k = ig[i], end = ig[i + 1]; k < end; k++)
		{
			pr[jg[k]] -= au[k] * val;
			al[k] = au[k] = 0;
		}
		for (int on_i = i + 1; on_i < n; on_i++) 
		{
			for (int k = ig[on_i], end = ig[on_i + 1]; k < end; k++)
			{
				if (jg[k] == i)
				{
					pr[on_i] -= al[k] * val;
					al[k] = au[k] = 0;
				}
			}
		}
	}

	// Получить правую часть
	block & rp(int i) {
		return pr[i];
	}

	// Получить диагональ
	block & diag(int i) {
		return di[i];
	}

	block element(int i, int j) {
		if (i == j) return di[i];
		bool eto_verh = i < j;
		if (eto_verh) std::swap(i, j);

		for (int k = ig[i], ig1 = ig[i + 1]; k < ig1; k++)
		{
			if (jg[k] == j) return eto_verh ? au[k] : al[k];
		}
		return 0;
	}

	block elementLU(int i, int j) {
		if (i == j) return gdi[i];
		bool eto_verh = i < j;
		if (eto_verh) std::swap(i, j);

		for (int k = ig[i], ig1 = ig[i + 1]; k < ig1; k++)
		{
			if (jg[k] == j) return eto_verh ? ggu[k] : ggl[k];
		}
		return 0;
	}

	//void test(block * x) {
	//	block * r = new block[n];
	//	mul_matrix(al, au, x, r);
	//	real diff_norm = 0;
	//	real f_norm = 0;
	//	real norm = 0;
	//	for (int i = 0; i < n; ++i) {
	//		block diff = abs(r[i] - pr[i]);
	//		printf("(%+.8e, %+.8e) | (%+.8e, %+.8e) | (%+.8e, %+.8e)\n", r[i].real(), r[i].imag(), pr[i].real(), pr[i].imag(), diff.real(), diff.imag());
	//		norm += pow(diff.real(), 2) + pow(diff.imag(), 2);
	//	}
	//	printf("nev = %+.8e\n", sqrt(norm));
	//}

	void print(FILE * fp) {
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j) {
				block e = element(i, j);
				fprintf(fp, "%+.8e ", e);
			}
			fprintf(fp, "| %+.8e\n", pr[i]);
		}
		//for (int i = 0; i < n; ++i) {
		//	printf("%f ", di[i]);
		//}printf("\n");
		//for (int i = 0; i < m; ++i) {
		//	printf("%f ", al[i]);
		//}
	}

	void printLU(FILE * fp) {
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j) {
				fprintf(fp, "%f  ", elementLU(i, j));
			}
			fprintf(fp, "\n");
		}
	}

	int size() {
		return n;
	}

	real sign(real a) {
		return a < 0 ? -1 : 1;
	}

	real nev(block * x) {
		block * r = new block[n];
		mul_matrix(al, au, x, r);
		real diff_norm = 0;
		real f_norm = 0;
		for (int i = 0; i < n; ++i) {
			f_norm += pow(std::abs(pr[i]), 2);
			diff_norm += pow(std::abs(r[i] - pr[i]), 2);
		}
		return sqrt(diff_norm) / sqrt(f_norm);
	}

	//real scalar(block * a, block * b) {
	//	real r = 0;
	//	for (int i = 0; i < n; ++i) {
	//		r += a[i].real() * b[i].real() + a[i].imag() * b[i].imag();
	//	}
	//	return r;
	//}

	real scalar(block * a, block * b) {
		real r = 0;
		for (int i = 0; i < n; ++i) {
			r += a[i]*b[i];
		}
		return r;
	}

	void mul_matrix(block * _al, block * _au, block * _x, block * res)
	{
		for (int i = 0; i < n; i++)
		{
			res[i] = di[i] * _x[i];
			for (int j = ig[i], end = ig[i + 1]; j < end; j++)
			{
				res[jg[j]] += _au[j] * _x[i];
				res[i] += _al[j] * _x[jg[j]];
			}
		}
	}

	void sol_low(block * _al, block * _f)
	{
		if (_al == NULL)
		{
			for (int i = 0; i < n; i++) _f[i] /= gdi[i];
		}
		else
		{
			for (int i = 0; i < n; i++)
			{
				for (int j = ig[i], end = ig[i + 1]; j < end; j++)
				{
					_f[i] -= _al[j] * _f[jg[j]];
				}
				_f[i] /= gdi[i];
			}
		}
	}

	void sol_up(block * _au, block * _f)
	{
		if (_au == NULL)
		{
			for (int i = n - 1; i >= 0; i--) _f[i] /= gdi[i];
		}
		else
		{
			for (int i = n - 1; i >= 0; i--)
			{
				_f[i] /= gdi[i];
				for (int j = ig[i], end = ig[i + 1]; j < end; j++)
				{
					_f[jg[j]] -= _au[j] * _f[i];
				}

			}
		}
	}

	void lu_fac()
	{
		ggu = new block[m];
		ggl = new block[m];
		gdi = new block[n];

		for (int i = 0; i < n; i++)
		{
			gdi[i] = di[i];
			for (int j = ig[i], end = ig[i + 1]; j < end; j++)
			{
				block sum_l = 0;
				block sum_u = 0;
				int k = jg[j];
				int on_k = ig[k];
				int on_i = ig[i];
				int end_k = ig[k + 1];
				while (on_k < end_k && on_i < j)
				{
					if (jg[on_i] == jg[on_k])
					{
						sum_l += ggl[on_i] * ggu[on_k];
						sum_u += ggu[on_i] * ggl[on_k];
						on_i++; on_k++;
					}
					else if (jg[on_i] < jg[on_k]) on_i++;
					else on_k++;
				}
				ggl[j] = (al[j] - sum_l) / gdi[k];
				ggu[j] = (au[j] - sum_u) / gdi[k];
				gdi[i] -= ggu[j] * ggl[j];
			}
			if (is_bad(gdi[i]))
			{
				printf("Impossible to make the LU(sq) decomposition of the matrix %lf.\n", gdi[i]);
				exit(1);
			}
			gdi[i] = sqrt(gdi[i]);

		}
		//print();
	}


	int los(block * x)
	{
		block * r = new block[n];
		block * z = new block[n];
		block * p = new block[n];
		block * tarr = new block[n];
		block * tezz = new block[n];
		real alpha, betta, pp;

		//    out_vec(cout, "al: ", al, k);
		//    out_vec(cout, "au: ", au, k);
		//    out_vec(cout, "di: ", di, n);
		//    out_vec(cout, "f: ", f, n);
		// Начальные значения
		int it = 0;
		mul_matrix(al, au, x, z);
		for (int i = 0; i < n; i++) r[i] = this->pr[i] - z[i];
		sol_low(ggl, r);
		for (int i = 0; i < n; i++) z[i] = r[i];
		sol_up(ggu, z);
		mul_matrix(al, au, z, p);
		sol_low(ggl, p);
		real nev = scalar(r, r);
		// Итерационный процесс
		while (nev > EPS && it < maxiter)
		{
			it++;
			pp = scalar(p, p);
			if (pp < EPS) { printf("pp = 0\n"); break; }
			alpha = scalar(p, r) / pp;
			for (int i = 0; i < n; i++)
			{
				x[i] += alpha * z[i];
				r[i] -= alpha * p[i];
				tarr[i] = r[i];
			}
			sol_up(ggu, tarr);
			mul_matrix(al, au, tarr, tezz);
			sol_low(ggl, tezz);
			for (int i = 0; i < n; i++) { tarr[i] = r[i]; }
			sol_up(ggu, tarr);
			betta = -scalar(p, tezz) / pp;
			for (int i = 0; i < n; i++)
			{
				z[i] = tarr[i] + betta * z[i];
				p[i] = tezz[i] + betta * p[i];
			}
			nev -= alpha * alpha * pp;
		}
		DEL(r); DEL(z); DEL(p); DEL(tarr); DEL(tezz);
		return it;
	}

	void mul_up(block * _au, block * _x)
	{
		if (_au == NULL)
		{
			for (int i = 0; i < n; i++) _x[i] = gdi[i] * _x[i];
		}
		else
		{
			for (int i = 0; i < n; i++)
			{
				_x[i] = gdi[i] * _x[i];
				for (int j = ig[i], end = ig[i + 1]; j < end; j++)
				{
					_x[jg[j]] += _au[j] * _x[i];
				}
			}
		}
	}
	template<class T>
	void alloc_arr(T * & arr, size_t size)
	{
		free_arr(arr);
		if (size > 0) arr = new T[size];
	}

	template<class T>
	void alloc_arr(T ** & arr, size_t si, size_t sj)
	{
		free_arr(arr, si);
		if (si > 0) arr = new T *[si];
		for (size_t i = 0; i < si; i++) if (sj > 0) arr[i] = new T[sj];
	}

	template<class T>
	void free_arr(T * &arr)
	{
		if (arr) delete[] arr;
		arr = NULL;
	}

	template<class T>
	void free_arr(T ** &arr, int size)
	{
		if (arr)
		{
			for (int i = 0; i < size; i++)
				delete[] arr[i];
			delete[] arr;
		}
		arr = NULL;
	}

	int gmres(block * x, size_t m)
	{
		block * r = new block[n];
		block * w = new block[n];
		block * tarr = new block[n];
		block ** v = NULL;
		real * d = new real[m + 1];
		real ** H = NULL;
		alloc_arr(H, m + 1, m);
		alloc_arr(v, n, m);

		//print();
		int it = 0;
		mul_matrix(al, au, x, r);
		//out_vec(cout, "r: ", r, n);
		mul_up(ggu, x);
		for (int i = 0; i < n; i++) r[i] = this->pr[i] - r[i];
		//out_vec(cout, "r: ", r, n);
		sol_low(ggl, r);
		//out_vec(cout, "r: ", r, n);

		real nev = sqrt(scalar(r, r));
		while (nev > EPS && it < maxiter)
		{
			for (int j = 0; j < n; j++) v[j][0] = r[j] / nev;

			for (int u = 0; u < m; u++)
			{
				for (int i = 0; i < n; i++) tarr[i] = v[i][u];
				sol_up(ggu, tarr);
				mul_matrix(al, au, tarr, w);
				sol_low(ggl, w);

				for (int l = 0; l <= u; l++)
				{
					H[l][u] = 0;
					//for (int i = 0; i < n; i++) H[l][u] += v[i][l].real() * w[i].real() + v[i][l].imag() * w[i].imag();
					for (int i = 0; i < n; i++) H[l][u] += v[i][l] * w[i];
					for (int i = 0; i < n; i++) w[i] -= H[l][u] * v[i][l];
				}
				for (int l = u + 1; l <= m; l++) H[l][u] = 0;

				real _H = H[u + 1][u] = sqrt(scalar(w, w));
				if (_H == 0)
				{
					m = u + 1;
					break;
				}
				if (u != m - 1)
				{
					for (int i = 0; i < n; i++) v[i][u + 1] = w[i] / _H;
				}
			}

			d[0] = nev;
			for (int i = 1; i <= m; i++) d[i] = 0;

			//        H[0][0] = 1, H[0][1] = 2, H[0][2] = 3;
			//        H[1][0] = 1, H[1][1] = 4, H[1][2] = 5;
			//        H[2][0] = 0, H[2][1] = 1, H[2][2] = 6;

			//        d[0] = 6, d[1] = 10, d[2] = 7;

			//        for (int i =0; i < m; i++) out_vec(cout, "H: ", H[i], m);

			//        for (int i = 1; i < m; i++)
			//        {
			//            real _H = H[i][i-1] / H[i-1][i-1];
			//            for (int j = i; j < m; j++)
			//                H[i][j] -= H[i-1][j] * _H;
			//            d[i] -= d[i-1] * _H;
			//        }


			for (int i = 0; i < m; i++)
			{
				real s = H[i + 1][i] / sqrt(pow(H[i][i], 2) + pow(H[i + 1][i], 2));
				real c = H[i][i] / sqrt(pow(H[i][i], 2) + pow(H[i + 1][i], 2));
				for (int j = i; j < m; j++)
				{
					real _H_i = H[i][j];
					real _H_i1 = H[i + 1][j];
					real _d_i = d[i];
					real _d_i1 = d[i + 1];

					H[i][j] = c * _H_i + s * _H_i1;
					H[i + 1][j] = -s * _H_i + c * _H_i1;
					d[i] = c * _d_i + s * _d_i1;
					d[i + 1] = -s * _d_i + c * _d_i1;
				}
			}

			//for (int i =0; i < m; i++) out_vec(cout, "H_up: ", H[i], m);

			//for (int i = 0; i < m; i++) out_vec(cout, "H: ", H[i], m);


			for (int i = m - 1; i >= 0; i--)
			{
				for (int j = m - 1; j > i; j--)
					d[i] -= d[j] * H[i][j];
				d[i] /= H[i][i];
			}

			//out_vec(cout, "d: ", d, m);

			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < m; j++)
					x[i] += d[j] * v[i][j];
			}

			for (int i = 0; i < n; i++) tarr[i] = x[i];
			sol_up(ggu, tarr);
			mul_matrix(al, au, tarr, r);
			for (int i = 0; i < n; i++) r[i] = this->pr[i] - r[i];
			sol_low(ggl, r);
			nev = sqrt(scalar(r, r));
			it++;
		}
		sol_up(ggu, x);

		DEL(r); DEL(w); DEL(d); DEL(tarr);
		free_arr(v, n); free_arr(H, m + 1);
		return it;
	}


	int bcg(block * x)
	{
		block * r = new block[n];
		block * z = new block[n];
		block * p = new block[n];
		block * s = new block[n];
		block * tarr = new block[n];
		block * tezz = new block[n];
		real alpha, betta, pr, _pr;

		int it = 0;
		mul_matrix(al, au, x, p);
		for (int i = 0; i < n; i++) r[i] = this->pr[i] - p[i];
		sol_low(ggl, r);
		for (int i = 0; i < n; i++)
		{
			z[i] = p[i] = s[i] = r[i];
		}
		sol_up(ggu, z);
		real nev = scalar(r, r);

		_pr = scalar(p, r);
		// Итерационный процесс
		while (nev > EPS && it < maxiter)
		{
			it++;
			mul_matrix(al, au, z, tarr);
			sol_low(ggl, tarr);
			alpha = _pr / scalar(s, tarr);
			for (int i = 0; i < n; i++)
			{
				x[i] += alpha * z[i];
				r[i] -= alpha * tarr[i];
				tezz[i] = s[i];
			}
			sol_up(ggl, tezz);
			mul_matrix(au, al, tezz, tarr);
			sol_low(ggu, tarr);
			for (int i = 0; i < n; i++)
			{
				p[i] -= alpha * tarr[i];
				tezz[i] = r[i];
			}
			pr = scalar(p, r);
			betta = pr / _pr;
			sol_up(ggu, tezz);
			for (int i = 0; i < n; i++)
			{
				z[i] = tezz[i] + betta * z[i];
				s[i] = p[i] + betta * s[i];
			}
			std::swap(pr, _pr);
			nev = scalar(r, r);
		}

		DEL(r); DEL(z); DEL(p); DEL(s); DEL(tarr); DEL(tezz);
		return it;
	}

	int solve(block * x) {
		for (int i = 0; i < n; ++i) {
			x[i] = 0;
		}

		
		lu_fac();
		//return los(x);
		return bcg(x);
		//return gmres(x, 4);
	}


	
};





