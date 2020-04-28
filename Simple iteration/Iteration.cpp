#include "Iteration.h"
using namespace std;

/// ��������� �������
void GenerateMatrix(vector<vector<float> >& matr, int sz) {
	assert(sz > 0);
	srand(time(0));
	// ������� ���������� ��� �������� ������� (����� ��������������)
	for (int i = 0; i < sz; ++i) {
		for (int j = 0; j < sz; ++j) {
			matr[i][j] = rand() % 5 - 4;
		}
	}
	// ��������� ������������ ��������
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
	// ��������� matr[0][0]
	for (int j = 1; j < sz; ++j) {
		sum += matr[0][j];
	}
	matr[0][0] = -sum + 1;
}

int Jacobi(int n, vector<vector<float> >& matr, vector<float>& f, vector<float>& x, float eps, int k_max, float w) {
	int k = 0;
	float max_dif = 2 * eps;

	// ��������� ������ x_0
	for (int i = 0; i < n; ++i) {
		x[i] = f[i] / matr[i][i];
	}

	while (max_dif >= eps && k <= k_max) {

		max_dif = 0;
		vector<float> prev(x);

		// ��������� ������ x(k+1)

		for (int i = 0; i < n; ++i) {
			x[i] = f[i];

			// ����� �� i (���������� �������� ���������� ��������)
			for (int j = 0; j < i; ++j) {
				x[i] -= matr[i][j] * prev[j];
			}

			// ����� ����� i (���������� �������� ���������� ��������)
			for (int j = i + 1; j < n; ++j) {
				x[i] -= matr[i][j] * prev[j];
			}

			// ����� �� ����������� a(i,i)
			x[i] /= matr[i][i];

			// ������������� max |x_i(k+1) - x_i(k)|
			max_dif = max(max_dif, fabs(x[i] - prev[i]));
		}

		// ����������� ����� ��������
		++k;
	}

	return k;
}

int Relaxation(int n, vector<vector<float> >& matr, vector<float>& f, vector<float>& x, float eps, int k_max, float w) {
	int k = 0;
	float max_dif = 2 * eps;

	// ��������� ������ x_0
	for (int i = 0; i < n; ++i) {
		x[i] = f[i] / matr[i][i];
	}

	while (max_dif >= eps && k <= k_max) {

		max_dif = 0;
		float last_val = 0;

		// ��������� ������ x(k+1)

		for (int i = 0; i < n; ++i) {
			last_val = x[i];
			x[i] = f[i];

			// ����� �� i (���������� �������� ������� ��������)
			for (int j = 0; j < i; ++j) {
				x[i] -= matr[i][j] * x[j];
			}

			// ����� ����� i (���������� �������� ���������� ��������)
			for (int j = i + 1; j < n; ++j) {
				x[i] -= matr[i][j] * x[j];
			}

			// ���������� ����������
			x[i] = x[i] * w / matr[i][i] + (1 - w) * last_val;

			// ������������� max |x_i(k+1) - x_i(k)|
			max_dif = max(max_dif, fabs(x[i] - last_val));
		}

		// ����������� ����� ��������
		++k;
	}

	return k;
}

/// ����� ������� �� �����
void PrintMatrix(int sz, vector<vector<float> >& matr, vector<float>& B) {
	for (int i = 0; i < sz; ++i) {
		for (int j = 0; j < sz; ++j) {
			cout << setw(7) << fixed << right << setprecision(2) << matr[i][j] << '\t';
		}
		cout << " |\t" << fixed << right << setprecision(2) << B[i] << '\n';
	}
}

/// ���������� ������ ����� ����������� �������
void Multiply(int sz, vector<vector<float> >& matr, vector<float>& X, vector<float>& B) {
	for (int i = 0; i < sz; ++i) {
		float temp = 0;
		for (int j = 0; j < sz; ++j) {
			temp += matr[i][j] * X[j];
		}
		B[i] = temp;
	}
}

void Solve(int n, vector<vector<float> >& A, vector<float>& X, vector<float>& B, float w, int (*method)(int , vector<vector<float> >& , vector<float>& , vector<float>&, float, int, float)) {
	vector<vector<float> > matr(A);
	vector<float> F(B);
	vector<float> ans(n);

	int k_max = 1000;
	float eps = 0.0001;

	// �������� ������ ����� (����� ���� ����������)
	int k = method(n, matr, F, ans, eps, k_max, w);
	if (k <= k_max) 
		cout << "��������� �������� ���������� �� " << k << "-� ��������.\n";
	else 
		cout << "��������� �������� �� ���������� �� " << k_max << " ��������.\n";

	cout << "\n\t������ ������������� �������:\n";
	cout << "X* = (";
	for (int i = 0; i < n - 1; ++i) {
		cout << setprecision(3) << ans[i] << ", ";
	}
	cout << ans[n - 1] << ")\n";
}

/// �������
int Lab5_Solution() {

	srand(time(NULL));
	int n = rand() % 3 + 10;
	const int m = 10;

	// ������� �������
	vector<vector<float> > A(n, vector<float>(n));
	GenerateMatrix(A, n);	

	// ������� ������ ������� � ������ ����� �������
	vector<float> X(n), B(n);
	for (int i = 0; i < n; ++i) {
		X[i] = m + i;
	}
	Multiply(n, A, X, B);

	cout << "\t���������:\n";
	cout << "n = " << n << '\n';
	cout << "m = " << m << '\n';
	cout << "\n\t��������������� �������: \n";
	PrintMatrix(n, A, B);

	cout << "\n\t������ ������� �������:\n";
	cout << "X = (";
	for (int i = 0; i < n - 1; ++i) {
		cout << X[i] << ", ";
	}
	cout << X[n - 1] << ")\n";
	
	////////////////////////////////////////

	cout << "\n\t\t\t����� �����\n\n";

	Solve(n, A, X, B, 1, &Jacobi);

	////////////////////////////////////////

	cout << "\n\t\t\t����� ����������\n\n";

	for (float w = 0.5; w <= 3; w += 0.5) {
		cout << "\n\t w = " << setprecision(1) << w << "\n\n";
		Solve(n, A, X, B, w, &Relaxation);
	}	

	return 0;
}