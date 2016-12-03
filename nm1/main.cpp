#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <algorithm>
#include <functional>
#include "Matrix.h"

#define SEPS 1.e-8
#define NUMDOF 4
using namespace std;

struct Point {
	double x;
	double f;
	double w;
};

// Читаем sum k^i, m - кол-во отрезков.
double sum_k(double k, double m)
{
	if (std::abs(k - 1.0) < SEPS) { return m; }
	return (std::pow(k, m) - 1.0) / (k - 1.0);
}

// Разбиваем область на отрезки.
void divide(std::vector<double> & v, double a, double b, double k, int n)
{
	v.resize(n);
	double h0 = (b - a) / sum_k(k, n - 1);
	double p = a;  //текущая точка
	double h = h0; // текущий отрезок
	for (int i = 0; i < n - 1; i++, p += h, h *= k) {
		v[i] = p;
	}
	v[n - 1] = b;
}

int main(int argc, char *argv[]) {

	const char * file = "in.txt";

	if (argc > 2) {
		file = argv[1];
	}

	ifstream fin(file);
	vector<Point> points;

	// Считываем точки
	int N;
	fin >> N;
	for (int i = 0; i < N; ++i) {
		Point p;
		fin >> p.x >> p.f >> p.w;
		points.push_back(p);
	}

	double a, b, k; //a начало b конец k коэф разрядки 
	int n;  // кол-во узлов
	fin >> a >> b >> n >> k;

	vector<double> grid;
	

	double beta;
	fin >> beta;


	// Начинаем решать

    // Распихиваем точки по отрезкам
	sort(points.begin(), points.end(), [](const Point &a, const Point &b) { return a.x < b.x; });
	//a = points[0].x - 0.1;
	//b = points[N - 1].x + 0.1;
	divide(grid, a, b, k, n);
	vector<set<int>> points_by_grid(n-1);
	int s = 0;
	for (int i = 0; i < N; ++i) {
		for (; s < n-1; ++s) {
			if (grid[s] <= points[i].x && points[i].x < grid[s + 1]) {
				break;
			}
		}
		if (s == n - 1) {
			cerr << "out of grid";
			return 1;
		}
		points_by_grid[s].insert(i);
	}


	auto getind = [](int elem, int * ind) {
		for (int i = 0; i < NUMDOF; ++i) {
			ind[i] = 2 * elem + i; 
		}
	};

	Matrix m;
	m.build(n - 1, NUMDOF);
	//m.build(n - 1, 2 * n, NUMDOF, getind);

	vector<function<double(double)>> phi = {
		[](double x) {return 1.0 - 3.0 * x * x + 2.0 * x * x * x; },
		[](double x) {return x - 2.0 * x * x + x * x * x; },
		[](double x) {return 3.0 * x * x - 2.0 * x * x * x; },
		[](double x) {return -x * x + x * x * x; }
	};

	for (int e = 0; e < n - 1; ++e) {
		int ind[4];
		double matr[NUMDOF * NUMDOF];
		double b[NUMDOF];
		double h = grid[e + 1] - grid[e];
		double x0 = grid[e];
		auto ksi = [phi, h, x0](int i, double x) { return (i % 2 == 0 ? 1.0 : h) * phi[i]((x - x0) / h); };

		getind(e, ind);
		for (int i = 0; i < NUMDOF; ++i) {
			for (int j = 0; j < NUMDOF; ++j) {
				matr[i*NUMDOF + j] = 0;
				for (auto & k : points_by_grid[e]) {
					auto & p = points[k];
					matr[i*NUMDOF + j] += p.w * ksi(i, p.x) * ksi(j, p.x);
				}
			}
			b[i] = 0;
			for (auto & k : points_by_grid[e]) {
				auto & p = points[k];
				b[i] += p.w * ksi(i, p.x) * p.f;
			}
		}
		m.add(1.0, matr, ind, NUMDOF);

		matr[0] = matr[10] = 12.0 / h / h / h;
		matr[1] = matr[3] = matr[4] = matr[12] = 6.0 / h / h;
		matr[2] = matr[8] = -matr[0];
		matr[5] = matr[15] = 4.0 / h;
		matr[6] = matr[9] = matr[11] = matr[14] = -matr[1];
		matr[7] = matr[13] = 2.0 / h;
 		
		m.add(beta, matr, ind, NUMDOF);
		m.add_f(1.0, b, ind, NUMDOF);
	}

	//m.LLt();
	//m.print();
	//for (int i = 0; i < m.size(); ++i) {
	//	for (int j = 0; j < m.size(); ++j) cout << m.element(i, j) << " ";  cout << "| " << m.rp(i);  cout << endl;
	//}

	double * q = new double[m.size()];
	m.solve(q);

	auto p = [&](double x) {
		int s = 0;
		for (; s < n - 1; ++s) {
			if (grid[s] <= x && x < grid[s + 1]) {
				break;
			}
		}
		if (s == n - 1) return 0.0;
		double h = grid[s + 1] - grid[s];
		double x0 = grid[s];
		auto ksi = [phi, h, x0](int i, double x) { return (i % 2 == 0 ? 1.0 : h) * phi[i]((x - x0) / h); };
		double res = 0;
		int ind[] = { 2 * s, 2 * s + 1, 2 * s + 2, 2 * s + 3 };
		for (int i = 0; i < NUMDOF; ++i) {
			res += q[ind[i]] * ksi(i, x);
		}
		return res;
	};

	//for (int i = 0; i < m.size(); ++i) {
	//	cout << q[i] << " ";
	//}
	//cout << endl;

	ofstream fout("out.txt");
	for (double x = a; x < b; x += 0.1) {
		fout << x << "\t\t\t" << p(x) << endl;
	}

	cout << "done" << endl;
	delete[] q;
	return 0;
}