#include "Sweep.h"

/// Генератор матрицы
void GenerateMatrix(vector<vector<float> >& matr, int N, int k, int m) {
	assert(N > 0);
	srand(time(0));
	// Генерация трехдиагональной матрицы
	for (int i = 1; i < N; ++i) {
		matr[i][i - 1] = -k;
		matr[i][i] = m + k + i - 1;
		matr[i][i + 1] = m + i - 1;
	}
	matr[0][0] = m;
	matr[0][1] = m - 1;
	matr[N][N] = m + k + N - 1;
	matr[N][N - 1] = -k;
}

/// Правая прогонка
bool RightSweep(int N, vector<vector<float> >& matr, vector<float>& f, vector<float>& alpha, vector<float>& beta) {
	alpha[1] = (-matr[0][1]) / matr[0][0];
	beta[1] = f[0] / matr[0][0];
	for (int i = 1; i <= N - 1; ++i) {
		alpha[i + 1] = (-matr[i][i + 1]) / (matr[i][i] - (-matr[i][i - 1]) * alpha[i]);
		beta[i + 1] = (f[i] + (-matr[i][i - 1]) * beta[i]) / (matr[i][i] - (-matr[i][i - 1]) * alpha[i]);
		matr[i - 1][i - 1] = 1;
		if (i > 1)
			matr[i - 1][i - 2] = 0;
		matr[i - 1][i] = -alpha[i];
		f[i - 1] = beta[i];
	}
	matr[N - 1][N - 1] = 1;
	matr[N - 1][N - 2] = 0;
	matr[N - 1][N] = -alpha[N];
	f[N - 1] = beta[N];
	beta[N + 1] = (f[N] + (-matr[N][N - 1]) * beta[N]) / (matr[N][N] - (-matr[N][N - 1]) * alpha[N]);
	matr[N][N - 1] = 0;
	matr[N][N] = 1;
	f[N] = beta[N + 1];
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
void CalculateAnswer(int N, vector<float>& alpha, vector<float>& beta, vector<float>& ans) {
	ans[N] = beta[N + 1];
	for (int i = N - 1; i >= 0; --i) {
		ans[i] = (alpha[i + 1]) * ans[i + 1] + beta[i + 1];
	}
}

/// Решение
int Solve() {
	cout << "\n\t\t\tМетод правой прогонки:\n\n";

	srand(time(NULL));
	int N = rand() % 3 + 9;

	vector<vector<float> > A(N + 1, vector<float>(N + 1, 0));

	const int m = 10;
	const int k = 2;

	// Генерируем матрицу
	GenerateMatrix(A, N, k, m);	

	// Создаем вектор решений
	vector<float> X(N + 1);
	for (int i = 0; i < N + 1; ++i) {
		X[i] = i + 1;
	}

	// Вычисление правой части расширенной матрицы
	vector<float> B(N + 1);
	for (int i = 0; i < N + 1; ++i) {
		float temp = 0;
		for (int j = 0; j < N + 1; ++j) {
			temp += A[i][j] * X[j];
		}
		B[i] = temp;
	}

	cout << "\tПараметры:\n";
	cout << "N + 1 = " << N + 1 << '\n';
	cout << "m = " << m << '\n';
	cout << "k = " << k << '\n';

	cout << "\n\tСгенерированная матрица: \n";
	for (int i = 0; i < N + 1; ++i) {
		for (int j = 0; j < N + 1; ++j) {
			cout << setw(5) << A[i][j] << '\t';
		}
		cout << " |\t" << B[i] << '\n';
	}

	vector<float> alpha(N + 2), beta(N + 2);

	RightSweep(N, A, B, alpha, beta);

	// Выводим матрицу
	cout << "\n\tМатрица после прямой прогонки: \n";
	for (int i = 0; i < N + 1; ++i) {
		for (int j = 0; j < N + 1; ++j) {
			cout << setw(7) << fixed << right << setprecision(2) << A[i][j] << '\t';
		}
		cout << " |\t" << B[i] << '\n';
	}

	vector<float> ans(N + 1);
	CalculateAnswer(N, alpha, beta, ans);

	cout << "\n\tВектор приближенного решения:\n";
	cout << "X* = (";
	for (int i = 0; i < N; ++i) {
		cout << ans[i] << ", ";
	}
	cout << ans[N] << ")\n";
	cout << "\n\tВектор точного решения:\n";
	cout << "X = (";
	for (int i = 0; i < N; ++i) {
		cout << X[i] << ", ";
	}
	cout << X[N] << ")\n";

	vector<float> tVec(N + 1);
	for (int i = 0; i < N + 1; ++i) {
		tVec[i] = X[i] - ans[i];
	}
	float relError = 100 * Norm(tVec) / Norm(X);

	cout << "\n\tОтносительная погрешность вычислений: " << setprecision(10) << relError << " %\n";
	
	return 0;
}