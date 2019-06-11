#include "Danilevski.h"

using namespace std;

/// Генератор матрицы
void GenerateMatrix(vector<vector<float> >& matr, int n) {

	assert(n > 0);
	srand(time(0));

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			matr[i][j] = rand() % 101 - 50;
		}
	}

}

/// Метод Данилевского
bool Danilevski(int n, vector<vector<float> >& matr, vector<vector<float> >& M) {

	for (int k = n - 1; k > 0; --k) {

		// вектор коэффициентов для матрицы M^-1
		vector<float> forM = matr[k];

		// ведущий элемент
		float lead = matr[k][k - 1];

		// если нерегулярный случай в строке		
		if (fabs(lead) < pow(10, -8)) {
			return false;
		}		

		// делим элементы (k-1)-го столбца на ведущий элемент
		for (int i = 0; i < n; ++i) {
			matr[i][k - 1] /= lead;
		}
		M[k][k - 1] = 1. / lead;

		// делаем нули в строке
		for (int j = 0; j < n; ++j) { 
			if (j == k - 1) continue;
			// сохраняем элемент матрицы M
			M[k][j] = -matr[k][j] / lead;
			// действуем, как в методе Гаусса
			float coeff = matr[k][j];
			for (int i = 0; i < n; ++i) {
				matr[i][j] -= matr[i][k - 1] * coeff;
			}					
		}

		// умножаем на матрицу M^(-1) слева (равносильно изменению (k-1)-й строки)
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

/// Вывод матрицы на экран
void PrintMatrix(int n, vector<vector<float> >& matr, ostream& os) {
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			os << setw(7) << fixed << right << setprecision(4) << matr[i][j] << '\t';
		}
		os << '\n';
	}
}

// След матрицы
float Trace(int n, vector<vector<float> >& matr) {
	float trace = 0;
	for (int i = 0; i < n; ++i) {
		trace += matr[i][i];
	}
	return trace;
}

/// Решение
int Lab6_Solution() {

	srand(time(NULL));
	int n = 4;
	
	ofstream fout("output.txt"); // для вывода матриц M в файл

	vector<vector<float> > A(n, vector<float>(n));

	cout << "\n\t\t\tМетод Данилевского\n";

	bool regenerated = false;

	while (true) {

		GenerateMatrix(A, n);		

		vector<vector<float> > matr(A);
		vector<vector<float> > M(n, vector<float>(n));

		// Вызываем метод Данилевского
		if (Danilevski(n, matr, M)) {
			
			// Если метод отработал, то выводим результаты
			// иначе - пересоздаем матрицу

			if (regenerated) {
				cout << "\tПо техническим причинам матрица A была перезаполнена\n";
			}

			cout << "\n\tИсходная матрица A: \n";
			fout << "\n\tИсходная матрица A: \n";
			PrintMatrix(n, A);	
			PrintMatrix(n, A, fout);

			cout << "\n\tСлед матрицы A: \n";
			fout << "\n\tСлед матрицы A: \n";
			cout << Trace(n, A) << '\n';
			fout << Trace(n, A) << '\n';

			cout << "\n\tПромежуточные матрицы M: \n";
			fout << "\n\tПромежуточные матрицы M: \n";
			for (int number = n - 1; number > 0; --number) { // по всем матрицам M
				cout << "\tM (" << number << "): \n\n";
				fout << "\tM (" << number << "): \n\n";
				for (int row = 0; row < n; ++row) { // строка
					for (int col = 0; col < n; ++col) { // столбец
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

			cout << "\n\tКаноническая форма Фробениуса: \n";
			fout << "\n\tКаноническая форма Фробениуса: \n";
			PrintMatrix(n, matr);
			PrintMatrix(n, matr, fout);

			cout << "\n\tКоэффициент p1: \n" << matr[0][0] << '\n';
			fout << "\n\tКоэффициент p1: \n" << matr[0][0] << '\n';

			break;
		}
		regenerated = true;
	}	

	return 0;
}