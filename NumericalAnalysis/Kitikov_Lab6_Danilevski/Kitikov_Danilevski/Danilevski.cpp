#include "Danilevski.h"

using namespace std;

/// ��������� �������
void GenerateMatrix(vector<vector<float> >& matr, int n) {

	assert(n > 0);
	srand(time(0));

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			matr[i][j] = rand() % 101 - 50;
		}
	}

}

/// ����� ������������
bool Danilevski(int n, vector<vector<float> >& matr, vector<vector<float> >& M) {

	for (int k = n - 1; k > 0; --k) {

		// ������ ������������� ��� ������� M^-1
		vector<float> forM = matr[k];

		// ������� �������
		float lead = matr[k][k - 1];

		// ���� ������������ ������ � ������		
		if (fabs(lead) < pow(10, -8)) {
			return false;
		}		

		// ����� �������� (k-1)-�� ������� �� ������� �������
		for (int i = 0; i < n; ++i) {
			matr[i][k - 1] /= lead;
		}
		M[k][k - 1] = 1. / lead;

		// ������ ���� � ������
		for (int j = 0; j < n; ++j) { 
			if (j == k - 1) continue;
			// ��������� ������� ������� M
			M[k][j] = -matr[k][j] / lead;
			// ���������, ��� � ������ ������
			float coeff = matr[k][j];
			for (int i = 0; i < n; ++i) {
				matr[i][j] -= matr[i][k - 1] * coeff;
			}					
		}

		// �������� �� ������� M^(-1) ����� (����������� ��������� (k-1)-� ������)
		for (int j = 0; j < n; ++j) {
			float sum = 0;
			for (int i = 0; i < n; ++i) {
				sum += matr[i][j] * forM[i];
			}
			matr[k - 1][j] = sum;
		}		
	}
	return true;
}

/// ����� ������� �� �����
void PrintMatrix(int n, vector<vector<float> >& matr, ostream& os) {
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			os << setw(7) << fixed << right << setprecision(4) << matr[i][j] << '\t';
		}
		os << '\n';
	}
}

// ���� �������
float Trace(int n, vector<vector<float> >& matr) {
	float trace = 0;
	for (int i = 0; i < n; ++i) {
		trace += matr[i][i];
	}
	return trace;
}

/// �������
int Lab6_Solution() {

	srand(time(NULL));
	int n = 4;
	
	ofstream fout("output.txt"); // ��� ������ ������ M � ����

	vector<vector<float> > A(n, vector<float>(n));

	cout << "\n\t\t\t����� ������������\n";

	bool regenerated = false;

	while (true) {

		GenerateMatrix(A, n);		

		vector<vector<float> > matr(A);
		vector<vector<float> > M(n, vector<float>(n));

		// �������� ����� ������������
		if (Danilevski(n, matr, M)) {
			
			// ���� ����� ���������, �� ������� ����������
			// ����� - ����������� �������

			if (regenerated) {
				cout << "\t�� ����������� �������� ������� A ���� �������������\n";
			}

			cout << "\n\t�������� ������� A: \n";
			fout << "\n\t�������� ������� A: \n";
			PrintMatrix(n, A);	
			PrintMatrix(n, A, fout);

			cout << "\n\t���� ������� A: \n";
			fout << "\n\t���� ������� A: \n";
			cout << Trace(n, A) << '\n';
			fout << Trace(n, A) << '\n';

			cout << "\n\t������������� ������� M: \n";
			fout << "\n\t������������� ������� M: \n";
			for (int number = n - 1; number > 0; --number) { // �� ���� �������� M
				cout << "\tM (" << number << "): \n\n";
				fout << "\tM (" << number << "): \n\n";
				for (int row = 0; row < n; ++row) { // ������
					for (int col = 0; col < n; ++col) { // �������
						if (number == row + 1) {
							cout << M[number][col] << '\t';
							fout << M[number][col] << '\t';
						}
						else {
							cout << ((row == col) ? "1\t" : "0\t");
							fout << ((row == col) ? "1\t" : "0\t");
						}
					}
					cout << '\n';
					fout << '\n';
				}
				cout << '\n';
				fout << '\n';
			}

			cout << "\n\t������������ ����� ����������: \n";
			fout << "\n\t������������ ����� ����������: \n";
			PrintMatrix(n, matr);
			PrintMatrix(n, matr, fout);

			cout << "\n\t����������� p1: \n" << matr[0][0] << '\n';
			fout << "\n\t����������� p1: \n" << matr[0][0] << '\n';

			break;
		}
		regenerated = true;
	}	

	return 0;
}