#include "Decomposition.h"

/// Генератор матрицы
void GenerateMatrix(vector<vector<float> >& matr, int sz, int k) {
	assert(sz > 0);
	srand(time(0));
	// Генерация симметрической матрицы
	for (int j = 0; j < sz; ++j) {
		for (int i = 0; i < j; ++i) {
			matr[i][j] = rand() % 5 - 4;
			matr[j][i] = matr[i][j];
		}
	}
	// Вычисляем диагональные элементы
	for (int i = 1; i < sz; ++i) {
		float sum = 0;
		for (int j = 0; j < sz; ++j) {
			if (j != i) {
				sum += matr[i][j];
			}
		}
		matr[i][i] = -sum;
	}
	float sum = 0;
	// Вычисляем matr[0][0]
	for (int j = 1; j < sz; ++j) {
		sum += matr[0][j];
	}
	matr[0][0] = -sum + pow(10, k);
}

/// Разложение матрицы
bool Decompose(int sz, vector<vector<float> >& matr) {
	vector<float> t(sz); // Массив для временного хранения вычисленных на предыдущем шаге элементов k-го столбца
	// Выполняем преобразования нижнего треугольника матрицы
	for (int k = 0; k < sz - 1; ++k) { // цикл по k
		for (int i = k + 1; i < sz; ++i) { // цикл по строкам ниже k-й
			t[i] = matr[i][k];
			matr[i][k] /= matr[k][k];
			for (int j = k + 1; j <= i; ++j) { // цикл по элементам строки (до i)
				matr[i][j] -= matr[i][k] * t[j];
			}
		}
	}
	// Проверяем, что определитель ненулевой
	for (int i = 0; i < sz; ++i) {
		if (matr[i][i] == 0)
			return false;
	}
	return true;
}

/// Вычисление кубической векторной нормы
float Norm(vector<float>& vec) {
	assert(vec.size() != 0);
	float ans = fabs(vec[0]);
	for (int i = 1; i < vec.size(); ++i) {
		ans = max(ans, fabs(vec[i]));
	}
	return ans;
}

/// Решение системы
void CalculateAnswer(int sz, vector<vector<float> >& matr, vector<float>& X, vector<float>& B, vector<float>& ans) {
	// Решаем уравнение Ly = B
	vector<float> Y(sz);
	Y[0] = B[0];
	for (int i = 1; i < sz; ++i) {
		float res = B[i];
		for (int j = 0; j < i; ++j) {
			res -= matr[i][j] * Y[j];
		}
		Y[i] = res;
	}
	// Решаем уравнение Dz = y
	vector<float> Z(sz);
	for (int i = 0; i < sz; ++i) {
		Z[i] = Y[i] / matr[i][i];
	}
	// Решаем уравнение Ltx = z
	ans[sz - 1] = Z[sz - 1];
	for (int i = sz - 2; i >= 0; --i) {
		float res = Z[i];
		for (int j = i + 1; j < sz; ++j) {
			res -= matr[j][i] * ans[j];
		}
		ans[i] = res;
	}
}

/// Решение
int Solve() {
	cout << "\n\t\t\t*LDLt-разложение матриц*\n\n";

	srand(time(NULL));
	int sz = rand() % 3 + 10;

	vector<vector<float> > A(sz, vector<float>(sz));

	// Генерируем матрицу для k = 0
	GenerateMatrix(A, sz, 0);
	const int m = 10;

	// Создаем вектор решений
	vector<float> X(sz), B(sz);
	for (int i = 0; i < sz; ++i) {
		X[i] = m + i;
	}

	// Вычисление правой части расширенной матрицы
	for (int i = 0; i < sz; ++i) {
		float temp = 0;
		for (int j = 0; j < sz; ++j) {
			temp += A[i][j] * X[j];
		}
		B[i] = temp;
	}

	cout << "\tПараметры:\n";
	cout << "n = " << sz << '\n';
	cout << "m = " << m << '\n';
	cout << "k = " << 0 << '\n';

	cout << "\n\tСгенерированная матрица: \n";
	for (int i = 0; i < sz; ++i) {
		for (int j = 0; j < sz; ++j) {
			cout << setw(5) << A[i][j] << '\t';
		}
		cout << " |\t" << B[i] << '\n';
	}

	vector<vector<float> > matr(A);

	cout << "\n\tРезультат разложения:\n";
	if (Decompose(sz, matr) == true) {
		// Выводим матрицу
		for (int i = 0; i < sz; ++i) {
			for (int j = 0; j <= i; ++j) {
				cout << setw(7) << fixed << right << setprecision(2) << matr[i][j] << '\t';
			}
			cout << '\n';
		}
		vector<float> ans(sz);
		CalculateAnswer(sz, matr, X, B, ans);

		cout << "\n\tВектор приближенного решения:\n";
		cout << "X* = (";
		for (int i = 0; i < sz - 1; ++i) {
			cout << ans[i] << ", ";
		}
		cout << ans[sz - 1] << ")\n";
		cout << "\n\tВектор точного решения:\n";
		cout << "X = (";
		for (int i = 0; i < sz - 1; ++i) {
			cout << X[i] << ", ";
		}
		cout << X[sz - 1] << ")\n";

		vector<float> tVec(sz);
		for (int i = 0; i < sz; ++i) {
			tVec[i] = X[i] - ans[i];
		}
		float relError = 100 * Norm(tVec) / Norm(X);

		cout << "\n\tОтносительная погрешность вычислений: " << setprecision(10) << relError << " %\n";
	}
	else {
		cout << "Матрица системы вырожденная. Разложение невозможно!\n";
	}		
	
	////////////////////////////////////////

	cout << "\n\tПараметры:\n";
	cout << "n = " << sz << '\n';
	cout << "m = " << m << '\n';
	cout << "k = " << 2 << '\n';

	// Изменяем матрицу для k = 2
	A[0][0] = A[0][0] - 1 + pow(10, -2);

	// Вычисление правой части расширенной матрицы
	for (int i = 0; i < sz; ++i) {
		float temp = 0;
		for (int j = 0; j < sz; ++j) {
			temp += A[i][j] * X[j];
		}
		B[i] = temp;
	}

	cout << "\n\tРезультат разложения:\n";
	if (Decompose(sz, A) == true) {
		// Выводим матрицу
		for (int i = 0; i < sz; ++i) {
			for (int j = 0; j <= i; ++j) {
				cout << setw(7) << fixed << right << setprecision(2) << A[i][j] << '\t';
			}
			cout << '\n';
		}
		vector<float> ans(sz);
		CalculateAnswer(sz, A, X, B, ans);

		cout << "\n\tВектор приближенного решения:\n";
		cout << "X* = (";
		for (int i = 0; i < sz - 1; ++i) {
			cout << ans[i] << ", ";
		}
		cout << ans[sz - 1] << ")\n";
		cout << "\n\tВектор точного решения:\n";
		cout << "X = (";
		for (int i = 0; i < sz - 1; ++i) {
			cout << X[i] << ", ";
		}
		cout << X[sz - 1] << ")\n";

		vector<float> tVec(sz);
		for (int i = 0; i < sz; ++i) {
			tVec[i] = X[i] - ans[i];
		}
		float relError = 100 * Norm(tVec) / Norm(X);

		cout << "\n\tОтносительная погрешность вычислений: " << setprecision(10) << relError << " %\n";
	}
	else {
		cout << "Матрица системы вырожденная. Разложение невозможно!\n";
	}
	
	return 0;
}