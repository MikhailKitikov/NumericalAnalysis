#pragma once
#include "Libs.h"
using namespace std;

/// сигнум
int sgn(float val) { return (val > 0 ? 1 : (val < 0 ? -1 : 0)); }

class SymmetricMatrix {

private:

	int n;
	vector<vector<float> > matr;

	Vector u, v;
	Vector u_29, u_30, v_30, v_31;
	Vector u_49, u_50, v_50, v_51;
	float lambda1_1, lambda1_2;
	float lambda2_1, lambda2_2, lambda2_3;
	float vnorm;
	int maxind;

public:

	/// конструктор
	SymmetricMatrix(int n) {
		assert(n > 0);
		this->n = n;
		matr.resize(n);
		for (int i = 0; i < n; ++i)
			matr[i].resize(n);
		Generate();

		u = Vector(n, 0);
		v.resize(n);
	}

	/// заполнение матрицы
	void Generate() {
		srand(time(0));
		// Генерация симметрической матрицы
		for (int j = 0; j < n; ++j) {
			for (int i = 0; i < j; ++i) {
				matr[i][j] = rand() % 5 - 4;
				matr[j][i] = matr[i][j];
			}
		}
		// Вычисляем диагональные элементы
		for (int i = 0; i < n; ++i) {
			float sum = 0;
			for (int j = 0; j < n; ++j) {
				if (j != i) {
					sum += matr[i][j];
				}
			}
			matr[i][i] = -sum;
			// Вычисляем matr[0][0]
			if (i == 0) ++matr[i][i];
		}
	}

	/// доступ к элементу
	float get(int i, int j) {
		assert(i >= 0 && i < n && j >= 0 && j < n);
		return matr[i][j];
	}

	/// вывод на экран
	void Print(ostream& os = cout) {
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j) {
				os << setw(3) << matr[i][j] << '\t';
			}
			os << '\n';
		}
	}

	/// умножение справа на вектор
	Vector& multiplyByVector(Vector& arg) {
		Vector* res = new Vector(n, 0);
		for (int i = 0; i < n; ++i) {
			float sum = 0;
			for (int j = 0; j < n; ++j) {
				sum += matr[i][j] * arg[j];
			}
			(*res)[i] = sum;
		}
		return *res;
	}

	/// Нахождение собственных значений
	void PowerMethod() {

		cout << "\n\n\t\tПервое max по модулю собственное значение\n\n\n";
		FirstEigenValue();

		cout << "\n\n\t\tВторое max по модулю собственное значение\n\n\n";
		SecondEigenValue();

	}

private:

	/// Нахождение 1-го max по модулю собств. значения
	void FirstEigenValue() {

		u[0] = 1; // u(k)
		v = this->multiplyByVector(u); // v(k+1)

		for (int k = 1; k <= 50; ++k) {

			u = v / v.cubicNorm();
			v = this->multiplyByVector(u);

			// lambda 1
			maxind = v.maxIndex();
			lambda1_1 = v[maxind] * sgn(u[maxind]);			

			// lambda 1 (для симметричных)
			lambda1_2 = ScalarProduct(v, u) / ScalarProduct(u, u);

			// сохраняем нужное для lambda 2
			if (k == 29) {
				u_29 = u;
				v_30 = v;
			}
			if (k == 30) {
				u_30 = u;
				v_31 = v;
			}
			if (k == 49) {
				u_49 = u;
				v_50 = v;
			}

			// вывод необходимой информации на последних итерациях
			if (k >= 46 && k <= 50) {
				cout << "\n\tIter " << k << endl;
				cout << "u(k) = " << u << endl;
				cout << "lambda_1_1 = " << setprecision(5) << lambda1_1 << endl;
				cout << "lambda_1_2 = " << setprecision(5) << lambda1_2 << endl;
			}
		}

		u_50 = u;
		v_51 = v;

		cout << "\nlambda_1_1 = " << lambda1_1 << endl;
		cout << "|| v(51) - lambda1_1 * u(50) || = " << (v - u * lambda1_1).cubicNorm() << endl;

		cout << "lambda_1_2 = " << lambda1_2 << endl;
		cout << "|| v(51) - lambda1_2 * u(50) || = " << (v - u * lambda1_2).cubicNorm() << endl;
	}

	/// Нахождение 2-го max по модулю собств. значения
	void SecondEigenValue() {

	// случай 1

		maxind = (v_30 - u_29 * lambda1_1).maxIndex();
		lambda2_1 = (v_31[maxind] * v_30.cubicNorm() - lambda1_1 * v_30[maxind]) / (v_30[maxind] - lambda1_1 * u_29[maxind]);
		cout << "\n1) lambda2 = " << lambda2_1 << endl;

		Vector eigen2 = v_31 - u_30 * lambda1_1;
		cout << "u = " << eigen2 << endl;		

		cout << "||A * u – lambda2 * u|| = " << setprecision(5) << (this->multiplyByVector(eigen2) - eigen2 * lambda2_1).cubicNorm() << endl;

	// случай 2

		maxind = (v_50 - u_49 * lambda1_1).maxIndex();
		lambda2_2 = (v_51[maxind] * v_50.cubicNorm() - lambda1_1 * v_50[maxind]) / (v_50[maxind] - lambda1_1 * u_49[maxind]);
		cout << "\n2) lambda2_2 = " << lambda2_2 << endl;

		eigen2 = v_51 - u_50 * lambda1_1;
		cout << "u = " << eigen2 << endl;

		cout << "||A * u – lambda2 * u|| = " << setprecision(5) << (this->multiplyByVector(eigen2) - eigen2 * lambda2_2).cubicNorm() << endl;

	// случай 3

		maxind = (v_50 - u_49 * lambda1_2).maxIndex();
		lambda2_3 = (v_51[maxind] * v_50.cubicNorm() - lambda1_2 * v_50[maxind]) / (v_50[maxind] - lambda1_2 * u_49[maxind]);
		cout << "\n3) lambda2 = " << setprecision(5) << lambda2_3 << endl;

		eigen2 = v_51 - u_50 * lambda1_2;
		cout << "u = " << eigen2 << endl;

		cout << "||A * u – lambda2 * u|| = " << setprecision(5) << (this->multiplyByVector(eigen2) - eigen2 * lambda2_3).cubicNorm() << endl;
	}

};