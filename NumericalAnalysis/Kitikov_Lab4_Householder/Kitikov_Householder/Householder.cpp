#include "Householder.h"
using namespace std;

/// ��������� �������
void GenerateMatrix(vector<vector<float> >& matr, int sz, int k) {
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
	matr[0][0] = -sum + pow(10, -k);
}

/// ����� ������ ��� ������ �������� ��������
void GaussWithoutSelection(int sz, vector<vector<float> >& matr, vector<float>& B) {
	for (int k = 0; k < sz - 1; ++k) { // ��������� k-� ������
		for (int i = k + 1; i < sz; ++i) { // �������� �� ���� ����� ���� k-� k-� � ������ �������������
			float coeff = matr[i][k] / matr[k][k];
			for (int j = k; j < sz; ++j) { // ���� �� ��������� ���������� ������
				matr[i][j] -= matr[k][j] * coeff;
			}
			B[i] -= B[k] * coeff;
		}
		if (k == 0) {
			cout << "\n\t������� ����� ������� ����:\n";
			PrintMatrix(sz, matr, B);
		}
	}
	cout << "\n\t�������� �������:\n";
	PrintMatrix(sz, matr, B);
}

/// ����� ���������
void Householder(int n, vector<vector<float> >& matr, vector<float>& B, vector<float>& ans) {

	vector<float> diag(n); // ��� �������� ������������ ���������
	for (int i = 0; i < n; ++i) {
		diag[i] = matr[i][i];
	}

	// ������� w_k ����� ������� �� ����� ������� ������������ ������� (������� ���������)

	for (int k = 0; k < n - 1; ++k) {

	/// ��������� w_k

		/// ���� �������� ��������

		float sgn = ((diag[k] == 0) ? 0 : diag[k] / fabs(diag[k]));

		/// ||s_k||

		float norm_s = 0;
		for (int j = k; j < n; ++j) {
			norm_s += matr[j][k] * matr[j][k];
		}
		norm_s = sqrt(norm_s);

		/// ||s_k - sgn(a[k][k]) * ||s_k|| * e_k||

		float norm_denom;
		if (sgn == 0)
			norm_denom = (diag[k] + norm_s) * (diag[k] + norm_s);
		else
			norm_denom = (diag[k] + sgn * norm_s) * (diag[k] + sgn * norm_s);

		for (int j = k + 1; j < n; ++j) {
			norm_denom += matr[j][k] * matr[j][k];
		}
		norm_denom = sqrt(norm_denom);

		/// w_k = (s_k + sgn(a[0][0]) * ||s_k|| * e_k) / ||s_k + sgn(a[0][0]) * ||s_k|| * e_k||

		if (sgn == 0)
			matr[k][k] = (diag[k] + norm_s) / norm_denom;
		else
			matr[k][k] = (diag[k] + sgn * norm_s) / norm_denom;

		for (int j = k + 1; j < n; ++j) {			
			matr[j][k] /= norm_denom;
		}

	/// ��������� �������� �������

		float scalar_prod;

		diag[k] = -sgn * norm_s;

		for (int j = k + 1; j < n; ++j) { /// ���� �� ��������

			/// scalar_prod = (matr[k:n][j], w_k)

			scalar_prod = 0;
			for (int i = k; i < n; ++i) {
				scalar_prod += matr[i][j] * matr[i][k];
			}

			/// matr[k:n][j] = matr[k:n][j] - 2 * (matr[k:n][j], w_k) * w_k

			for (int i = k; i < n; ++i) {
				matr[i][j] -= 2 * scalar_prod * matr[i][k];
				if (i == j) {
					diag[i] = matr[i][i]; // ��������� ������������ �������
				}
			}
		}

	/// ��������� ������ �����

		scalar_prod = 0;
		for (int i = k; i < n; ++i) {
			scalar_prod += B[i] * matr[i][k];
		}

		/// ������ �����: B = B - 2 * (B, w_k) * w_k

		for (int i = k; i < n; ++i) {
			B[i] -= 2 * scalar_prod * matr[i][k];
		}

		if (k == 0) {
			cout << "\n\t������� ����� ������� ����:\n";
			PrintMatrix(n, matr, B);
			cout << "\n\t������������ ��������:\n";
			for (int i = 0; i < n; ++i) {
				cout << setprecision(2) << diag[i] << '\t';
			}
			cout << endl;
		}
	}

	/// ��������� �������

	ans[n - 1] = B[n - 1] / diag[n - 1];
	for (int i = n - 2; i >= 0; --i) {
		float res = B[i];
		for (int j = i + 1; j < n; ++j) {
			res -= matr[i][j] * ans[j];
		}
		res /= diag[i];
		ans[i] = res;
	}

	cout << "\n\t�������� �������:\n";
	PrintMatrix(n, matr, B);
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

/// ���������� ���������� ��������� �����
float CubicNorm(vector<float>& vec) {
	assert(vec.size() != 0);
	float ans = fabs(vec[0]);
	for (int i = 1; i < vec.size(); ++i) {
		ans = max(ans, fabs(vec[i]));
	}
	return ans;
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

/// �������� ��� ������ ������
void CalculateAnswer(int sz, vector<vector<float> >& matr, vector<float>& X, vector<float>& B, vector<float>& ans) {
	ans[sz - 1] = B[sz - 1] / matr[sz - 1][sz - 1];
	for (int i = sz - 2; i >= 0; --i) {
		float res = B[i];
		for (int j = i + 1; j < sz; ++j) {
			res -= matr[i][j] * ans[j];
		}
		res /= matr[i][i];
		ans[i] = res;
	}
}

/// �������
int Solve() {
	cout << "\n\t\t\t����� ������ ��� ������ �������� ��������\n\n";

	srand(time(NULL));
	int sz = rand() % 3 + 10;

	vector<vector<float> > A(sz, vector<float>(sz));

	// ���������� ������� ��� k = 2
	const int k = 2;
	GenerateMatrix(A, sz, k);
	const int m = 10;

	vector<float> X(sz), B(sz);
	for (int i = 0; i < sz; ++i) {
		X[i] = m + i;
	}

	Multiply(sz, A, X, B);

	cout << "\t���������:\n";
	cout << "n = " << sz << '\n';
	cout << "m = " << m << '\n';
	cout << "k = " << k << '\n';

	cout << "\n\t��������������� �������: \n";
	for (int i = 0; i < sz; ++i) {
		for (int j = 0; j < sz; ++j) {
			cout << setw(5) << A[i][j] << '\t';
		}
		cout << " |\t" << B[i] << '\n';
	}

	vector<vector<float> > matr(A);
	// � ������ ������ �������� � �������� matr - ������ ������� A
	// � � ������ ��������� - � ����� A

	GaussWithoutSelection(sz, matr, B);

	vector<float> ans(sz);
	CalculateAnswer(sz, matr, X, B, ans);

	cout << "\n\t������ ������������� �������:\n";
	cout << "X* = (";
	for (int i = 0; i < sz - 1; ++i) {
		cout << ans[i] << ", ";
	}
	cout << ans[sz - 1] << ")\n";
	cout << "\n\t������ ������� �������:\n";
	cout << "X = (";
	for (int i = 0; i < sz - 1; ++i) {
		cout << X[i] << ", ";
	}
	cout << X[sz - 1] << ")\n";

	vector<float> tVec(sz);
	for (int i = 0; i < sz; ++i) {
		tVec[i] = X[i] - ans[i];
	}
	float relError = 100 * CubicNorm(tVec) / CubicNorm(X);

	cout << "\n\t������������� ����������� ����������: " << setprecision(10) << relError << " %\n";
	
	////////////////////////////////////////

	cout << "\n\t\t\t����� ���������\n\n";

	for (int i = 0; i < sz; ++i) {
		X[i] = m + i;
	}
	Multiply(sz, A, X, B);

	Householder(sz, A, B, ans);

	cout << "\n\t������ ������������� �������:\n";
	cout << "X* = (";
	for (int i = 0; i < sz - 1; ++i) {
		cout << ans[i] << ", ";
	}
	cout << ans[sz - 1] << ")\n";
	cout << "\n\t������ ������� �������:\n";
	cout << "X = (";
	for (int i = 0; i < sz - 1; ++i) {
		cout << X[i] << ", ";
	}
	cout << X[sz - 1] << ")\n";

	for (int i = 0; i < sz; ++i) {
		tVec[i] = X[i] - ans[i];
	}
	relError = 100 * CubicNorm(tVec) / CubicNorm(X);

	cout << "\n\t������������� ����������� ����������: " << setprecision(10) << relError << " %\n";

	return 0;
}