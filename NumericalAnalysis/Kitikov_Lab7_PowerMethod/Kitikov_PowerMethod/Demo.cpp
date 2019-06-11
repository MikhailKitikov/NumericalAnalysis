#include "Libs.h"
using namespace std;

void Solve() {

	cout << "\n\t\t\t*Итерационный степенной метод*\n\n";

	srand(time(NULL));
	int n = rand() % 3 + 10;	
	cout << "\nn = " << n << endl;
	SymmetricMatrix A(n);
	cout << "\n\tСгенерированная матрица: \n";
	A.Print();

	A.PowerMethod();
}

int main() {

	setlocale(LC_ALL, ".1251");

	Solve();

	system("pause");
	return 0;
}
