#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <algorithm>
#include <functional>
#include <random>
#include <chrono>
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

void generate(const char * file, function<double(double)> f, double a, double b, int n, double sig = 0.1, double sigout = 1.0, double out = 0.3, double w = 1.0) {
	vector<double> grid;
	divide(grid, a, b, 1.0, n);
	ofstream fout(file);
	ofstream fclean("clean.txt");
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	default_random_engine generator(seed);
	normal_distribution<double> normal_dist(0, sig);
	normal_distribution<double> out_dist(0, sigout);
	binomial_distribution<int> binom_dist(1, out);
	auto normal = bind(normal_dist, generator);
	auto binom  = bind(binom_dist, generator);
	auto outlier = bind(out_dist, generator);
	fout << grid.size() << endl;
	fclean << grid.size() << endl;
	for (int i = 0; i < grid.size(); ++i) {
		double x = grid[i];
		double err = binom() ? outlier() : normal();
		fout << x << "\t" << f(x) + err << "\t" << w << endl;
		fclean << x << "\t" << f(x) << "\t" << w << endl;
	}
}


int main(int argc, char *argv[]) {

	const char * file = "in.txt";
	const char * filepoints = "points.txt";
	//const char * filepoints = "clean.txt";

	if (argc > 2) {
		file = argv[1];
	}

	//generate(filepoints, [](double x) {return sin(x); }, -3.5, 3.5, 15);

	ifstream fin(filepoints);
	vector<Point> points;

	// Считываем точки
	int N;
	fin >> N;
	for (int i = 0; i < N; ++i) {
		Point p;
		fin >> p.x >> p.f >> p.w;
		points.push_back(p);
	}

	fin.close();
	fin.open(file);

	double a, b, k; //a начало b конец k коэф разрядки 
	int n;  // кол-во узлов
	fin >> n >> k;

	vector<double> grid;
	

	double beta;
	fin >> beta;

	double r;  // заданное кол-во раз
	fin >> r;

	// Начинаем решать

    // Распихиваем точки по отрезкам
	sort(points.begin(), points.end(), [](const Point &a, const Point &b) { return a.x < b.x; });
	double h = points.back().x - points.front().x;
	double out = 0.05;
	a = points.front().x - h * out;
	b = points.back().x + h * out;
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
	vector<double> q(m.size());
	vector<function<double(double)>> phi = {
		[](double x) {return 1.0 - 3.0 * x * x + 2.0 * x * x * x; },
		[](double x) {return x - 2.0 * x * x + x * x * x; },
		[](double x) {return 3.0 * x * x - 2.0 * x * x * x; },
		[](double x) {return -x * x + x * x * x; }
	};

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

	vector<double> err(N);
	for (int it = 0; it< 200; ++it) { // цикл по фильтрации
		m.clean();
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
		m.solve(q);

		double mean = 0;
		double max = 0;
		

		for (int i = 0; i < N; ++i) {
			auto & point = points[i];
			err[i] = fabs(p(point.x) - point.f);
			if (err[i] > max) { max = err[i]; }
			mean += err[i];
		}
		mean /= N;

		cout << "it = " << it << "\tmean = " << mean << "\tmax = " << max << "\tratio = " << max / mean << endl;

		bool changed = false;
		for (int i = 0; i < N; ++i) {
			auto & point = points[i];
			if (err[i] > mean * r) {
				point.w *= 1.05;
				cout << "err[" << i << "] = " << err[i] << " > mean = " << mean << " in " << err[i] / mean << " times. w[" << i << "] set to " << point.w << endl;
				changed = true;
			}
		}
		if (!changed) {
			break;
		}
	}



	ofstream fout("out.txt");
	for (double x = a; x <= b; x += 0.1) {
		fout << x << "\t" << p(x) << endl;
	}

	cout << "done" << endl;
	return 0;
}