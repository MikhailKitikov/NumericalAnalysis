#include "Gauss.h"

/// ��������� ������� ��� ������� �1
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
	matr[0][0] = -sum + pow(10, k);
}

/// ����� ������ ��� ������ �������� ��������
void GaussWithoutSelection(int sz, vector<vector<float> >& matr, vector<float>& B) {
	// ������ 1-� ���
	for (int i = 1; i < sz; ++i) {
		if (matr[0][0] == 0) { throw "Single solution can't be found.\n"; }
		float coeff = matr[i][0] / matr[0][0];
		for (int j = 0; j < sz; ++j) {
			matr[i][j] -= matr[0][j] * coeff;
		}
		B[i] -= B[0] * coeff;
	}
	cout << "\n\t������� ����� ������� ����:\n";
	for (int i = 0; i < sz; ++i) {
		for (int j = 0; j < sz; ++j) {
			cout << setw(7) << fixed << right << setprecision(2) << matr[i][j] << '\t';
		}
		cout << " |\t" << fixed << right << setprecision(3) << B[i] << '\n';
	}
	// ��������� ����
	for (int k = 1; k < sz - 1; ++k) { // ���� �� k
		for (int i = k + 1; i < sz; ++i) { // ���� �� ������� ���� k-�
			if (matr[k][k] == 0) { throw "Single solution can't be found.\n"; }
			float coeff = matr[i][k] / matr[k][k];
			for (int j = k; j < sz; ++j) { // ���� �� ��������� ������
				matr[i][j] -= matr[k][j] * coeff;
			}
			B[i] -= B[k] * coeff;
		}
	}
	cout << "\n\t�������� �������:\n";
	for (int i = 0; i < sz; ++i) {
		for (int j = 0; j < sz; ++j) {
			cout << setw(7) << fixed << right << setprecision(2) << matr[i][j] << '\t';
		}
		cout << " |\t" << fixed << right << setprecision(3) << B[i] << '\n';
	}
}

/// ����� ������ � ������� �������� �������� �� �������
void GaussWithSelection(int sz, vector<vector<float> >& matr, vector<float>& B) {
	int maxind = 0;
	for (int i = 1; i < sz; ++i) { // ���� max � 1 �������
		if (fabs(matr[i][0]) > fabs(matr[maxind][0]))
			maxind = i;
	}
	cout << "\n\t�� 1-� ���� ������� ������� ���� matr[" << maxind << "][0] = " << matr[maxind][0] << '\n';
	if (maxind != 0) { // ������ ������� ������ (���� ����)
		swap(matr[0], matr[maxind]);
		swap(B[0], B[maxind]);
	}		
	for (int i = 1; i < sz; ++i) { // ������� ��� ������ ������
		if (matr[0][0] == 0) { throw "Single solution can't be found.\n"; }
		float coeff = matr[i][0] / matr[0][0];
		for (int j = 0; j < sz; ++j) {
			matr[i][j] -= matr[0][j] * coeff;
		}
		B[i] -= B[0] * coeff;
	}
	cout << "\n\t������� ����� ������� ����:\n";
	for (int i = 0; i < sz; ++i) {
		for (int j = 0; j < sz; ++j) {
			cout << setw(7) << fixed << right << setprecision(2) << matr[i][j] << '\t';
		}
		cout << " |\t" << fixed << right << setprecision(3) << B[i] << '\n';
	}
	for (int k = 1; k < sz - 1; ++k) { // ������ ��������� ����
		maxind = k;
		for (int i = k; i < sz; ++i) {
			if (fabs(matr[i][k]) > fabs(matr[maxind][k]))
				maxind = i;
		}
		if (maxind != k) {
			swap(matr[k], matr[maxind]);
			swap(B[k], B[maxind]);
		}			
		for (int i = k + 1; i < sz; ++i) {
			if (matr[k][k] == 0) { throw "Single solution can't be found.\n"; }
			float coeff = matr[i][k] / matr[k][k];
			for (int j = k; j < sz; ++j) {
				matr[i][j] -= matr[k][j] * coeff;
			}
			B[i] -= B[k] * coeff;
		}
	}
	cout << "\n\t�������� �������:\n";
	for (int i = 0; i < sz; ++i) {
		for (int j = 0; j < sz; ++j) {
			cout << setw(7) << fixed << right << setprecision(2) << matr[i][j] << '\t';
		}
		cout << " |\t" << fixed << right << setprecision(3) << B[i] << '\n';
	}
}

/// ���������� ���������� ��������� �����
float Norm(vector<float>& vec) {
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

/// ������� �1
int Task1() {
	cout << "\n\t\t\t\t\t������� �1\n\n";
	cout << "\n\t\t\t*����� ������ ��� ������ �������� ��������*\n\n";

	srand(time(NULL));
	int sz = rand() % 4 + 12;

	vector<vector<float> > A(sz, vector<float>(sz));

	// ���������� ������� ��� k = 0
	GenerateMatrix(A, sz, 0);
	const int m = 10;

	vector<float> X(sz), B(sz);
	for (int i = 0; i < sz; ++i) {
		X[i] = m + i;
	}

	Multiply(sz, A, X, B);

	cout << "\t���������:\n";
	cout << "n = " << sz << '\n';
	cout << "m = " << m << '\n';
	cout << "k = " << 0 << '\n';

	cout << "\n\t��������������� �������: \n";
	for (int i = 0; i < sz; ++i) {
		for (int j = 0; j < sz; ++j) {
			cout << setw(5) << A[i][j] << '\t';
		}
		cout << " |\t" << B[i] << '\n';
	}

	vector<vector<float> > matr(A);

	try {
		GaussWithoutSelection(sz, matr, B);
	}
	catch (const char* s) {
		cout << s;
		return 1;
	}

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
	float relError = 100 * Norm(tVec) / Norm(X);

	cout << "\n\t������������� ����������� ����������: " << setprecision(10) << relError << " %\n";
	
	////////////////////////////////////////

	cout << "\n\t���������:\n";
	cout << "n = " << sz << '\n';
	cout << "m = " << m << '\n';
	cout << "k = " << 2 << '\n';

	// �������� ������� ��� k = 2
	A[0][0] = A[0][0] - 1 + pow(10, -2);

	Multiply(sz, A, X, B);

	GaussWithoutSelection(sz, A, B);

	CalculateAnswer(sz, A, X, B, ans);

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
	relError = 100 * Norm(tVec) / Norm(X);

	cout << "\n\t������������� ����������� ����������: " << setprecision(10) << relError << " %\n";

	return 0;
}

////////////////////

/// ������� �2
int Task2() {
	cout << "\n\t\t\t\t\t\t������� �2\n\n";

	srand(time(NULL));
	int sz = rand() % 6 + 15;

	vector<vector<float> > A(sz, vector<float>(sz));

	// ���������� �������
	for (int i = 0; i < sz; ++i) {
		for (int j = 0; j < sz; ++j) {
			A[i][j] = rand() % 201 - 100;
		}
	}
	
	const int m = 10;

	vector<float> X(sz), B(sz);
	for (int i = 0; i < sz; ++i) {
		X[i] = m + i;
	}

	Multiply(sz, A, X, B);

	cout << "\t���������:\n";
	cout << "n = " << sz << '\n';
	cout << "m = " << m << '\n';

	cout << "\n\t��������������� �������: \n";
	for (int i = 0; i < sz; ++i) {
		for (int j = 0; j < sz; ++j) {
			cout << setw(5) << setprecision(0) << A[i][j] << '\t';
		}
		cout << " |\t" << B[i] << '\n';
	}

	cout << "\n\t\t\t*����� ������ ��� ������ �������� ��������*\n\n";

	vector<vector<float> > matr(A);

	try {
		GaussWithoutSelection(sz, matr, B);
	}
	catch (const char* s) {
		cout << s;
		return 1;
	}	

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
	float relError = 100 * Norm(tVec) / Norm(X);

	cout << "\n\t������������� ����������� ����������: " << setprecision(10) << relError << " %\n";

	////////////////////////

	cout << "\n\t\t\t*����� ������ � ������� �������� �������� �� �������*\n\n";

	Multiply(sz, A, X, B);

	try {
		GaussWithSelection(sz, A, B);
	}
	catch (const char* s) {
		cout << s;
		return 1;
	}	

	CalculateAnswer(sz, A, X, B, ans);

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
	relError = 100 * Norm(tVec) / Norm(X);

	cout << "\n\t������������� ����������� ����������: " << setprecision(10) << relError << " %\n";

	return 0;
}